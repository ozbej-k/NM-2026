module dn1
using LinearAlgebra, Graphs

export RedkaMatrika, sor, krožna_lestev, matrika, desne_strani, vloži!

"""
RedkaMatrika(V, I): Sparse square matrix definition where:
- V[i][j] stores non-zero elements from row i
- I[i][j] stores the column index for each non-zero element in V
"""
mutable struct RedkaMatrika{T}
    V::Vector{Vector{T}}
    I::Vector{Vector{Int}}
end

"""
Returns size of matrix.
"""
function Base.size(A::RedkaMatrika)
    # Return the number of rows, since the matrix is square, that is also the number of columns
    return (length(A.V), length(A.V))
end

Base.size(A::RedkaMatrika, d::Int) = size(A)[d]

"""
Returns value of element at index (i, j).
"""
function Base.getindex(A::RedkaMatrika{T}, i::Int, j::Int) where T
    @boundscheck 1 <= i <= length(A.V) || throw(BoundsError(A, (i, j)))
    @boundscheck 1 <= j <= length(A.V) || throw(BoundsError(A, (i, j)))

    # Search for the correct column
    for _j in eachindex(A.I[i])
        if A.I[i][_j] == j
            return A.V[i][_j]
        end
    end

    # No element found, value is zero
    return zero(T)
end # Base.getindex

"""
Sets value of element at index (i, j).
"""
function Base.setindex!(A::RedkaMatrika{T}, value, i::Int, j::Int) where T
    @boundscheck 1 <= i <= length(A.V) || throw(BoundsError(A, (i, j)))
    @boundscheck 1 <= j <= length(A.V) || throw(BoundsError(A, (i, j)))

    # Search for existing element
    for k in eachindex(A.I[i])
        if A.I[i][k] == j
            if value == zero(T)
                # Remove when setting non-zero to zero
                deleteat!(A.I[i], k)
                deleteat!(A.V[i], k)
            else
                A.V[i][k] = convert(T, value)
            end
            return A
        end
    end

    # Element does not exist, add new element
    if value != zero(T)
        push!(A.I[i], j)
        push!(A.V[i], convert(T, value))
    end
    return A
end # Base.setindex!

"""
Multiplying sparse matrix with vector (from right)
"""
function Base.:*(A::RedkaMatrika{T}, v::AbstractVector) where T
    n, _ = size(A)
    length(v) == n || throw(DimensionMismatch("Matrix has $n columns but vector has length $(length(v))"))

    y = zeros(T, n)
    for i in 1:n
        s = zero(T)
        for k in eachindex(A.V[i])
            s += A.V[i][k] * v[A.I[i][k]]
        end
        y[i] = s
    end
    return y
end # Base.*

"""
SOR iteration for sparse matrices, solves system Ax = b.
"""
function sor(A::RedkaMatrika{T}, b::AbstractVector{T}, x0::AbstractVector{T}, omega, tol=1e-10) where T
    n, _ = size(A)
    length(b) == n || throw(DimensionMismatch("Matrix is size ($n, $n) but b has length $(length(b))"))
    length(x0) == n || throw(DimensionMismatch("Matrix is size ($n, $n) but x0 has length $(length(x0))"))

    # SOR iteration
    x = copy(x0)
    it = 0
    while true
        it += 1
        for i in 1:n
            diag = zero(T)
            s = zero(T)

            for k in eachindex(A.V[i])
                j = A.I[i][k]
                aij = A.V[i][k]

                if j == i
                    diag = aij
                else
                    s += aij * x[j]
                end
            end
            x[i] = (1 - omega) * x[i] + omega * (b[i] - s) / diag
        end

        # Check convergence
        if norm(A * x - b, Inf) < tol
            return x, it
        end
    end
end # sor


"""
Function that creates a circular ladder graph.
"""
function krožna_lestev(n)
    G = SimpleGraph(2 * n)
    # prvi cikel
    for i = 1:n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, 1, n)
    # drugi cikel
    for i = n+1:2n-1
        add_edge!(G, i, i + 1)
    end
    add_edge!(G, n + 1, 2n)
    # povezave med obema cikloma
    for i = 1:n
        add_edge!(G, i, i + n)
    end
    return G
end # krožna_lestev

"""
Function that creates the matrix of the system for force equilibrium in a graph.
"""
function matrika(G::AbstractGraph, sprem)
    # preslikava med vozlišči in indeksi v matriki
    v_to_i = Dict([sprem[i] => i for i in eachindex(sprem)])
    m = length(sprem)
    
    A = RedkaMatrika([Float64[] for _ in 1:m], [Int[] for _ in 1:m])
    for i = 1:m
        vertex = sprem[i]
        sosedi = neighbors(G, vertex)
        for vertex2 in sosedi
            if haskey(v_to_i, vertex2)
                j = v_to_i[vertex2]
                A[i, j] = 1
            end
        end
        A[i, i] = -length(sosedi)
    end

    return A
end # matrika

"""
Function that calculates the right-hand sides of the system for force equilibrium in a graph 
based on the coordinates of the fixed vertices.
"""
function desne_strani(G::AbstractGraph, sprem, koordinate)
    set = Set(sprem)
    m = length(sprem)
    b = zeros(m)
    for i = 1:m
        v = sprem[i]
        for v2 in neighbors(G, v)
            if !(v2 in set) # dodamo le točke, ki so fiksirane
                b[i] -= koordinate[v2]
            end
        end
    end
    return b
end

"""
Function that finds the coordinates of a graph's embedding using the physical method.
"""
function vloži!(G::AbstractGraph, fix, točke, omega; tol=10e-10) 
    sprem = setdiff(vertices(G), fix)
    dim, _ = size(točke)
    A = matrika(G, sprem)
    neg_a = RedkaMatrika(-1.0 * A.V, A.I)

    n, _ = size(A)
    x0 = zeros(n)
    it_sum = 0
    for k = 1:dim
        b = desne_strani(G, sprem, točke[k, :])
        x, it = sor(neg_a, -b, x0, omega, tol) # matrika A je negativno definitna
        točke[k, sprem] = x
        it_sum += it
    end

    return it_sum
end # vloži!

end # module dn1