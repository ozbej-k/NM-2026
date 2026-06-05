module dn2

using LinearAlgebra

export force, my_gausslegendre

"""
Data type for closed interval [min, max]
"""
struct Interval{T}
    min::T
    max::T
end

"""
Data type for definite integral
"""
struct Integral{T}
    f
    interval::Interval{T}
end

"""
Length of interval
"""
length(int::Interval) = int.max - int.min

"""
Data type for the Gauss-Legendre quadrature with nodes x and weights u on given inverval
"""
struct Quadrature{T}
    x::Vector{T}
    u::Vector{T}
    interval::Interval{T}
end

"""
A function that calculates the integral using the Gauss-Legendre quadrature
"""
function integrate(i::Integral, k::Quadrature)
    scale = length(i.interval) / length(k.interval)
    map(x) = scale * (x - k.interval.min) + i.interval.min
    I = sum(u * i.f(map(x)) for (u, x) in zip(k.u, k.x))
    return scale * I
end

"""
A function that calculates the Gauss-Legendre quadrature nodes and weights
"""
function my_gausslegendre(n)
    diagonal = zeros(n)
    subdiagonal = [i / sqrt(4i^2 - 1) for i in 1:(n-1)]
    J = SymTridiagonal(diagonal, subdiagonal)
    x, V = eigen(J)
    w = 2.0 .* (V[1, :].^2)
    return x, w
end

"""
A function that returns Gauss-Legendre quadrature with n nodes
"""
function glkvad(n)
    x, u = my_gausslegendre(n)
    return Quadrature(x, u, Interval(-1.0, 1.0))
end


"""
A function which calculates the gravitational force between two parallel-placed unit homogeneous cubes at a distance of 1 
using the Gauss-Legendre quadrature with n nodes.
"""
function force(n)
    kv = glkvad(n)
    I01 = Interval(0.0, 1.0)
    I23 = Interval(2.0, 3.0)

    I_x1 = Integral(x1 -> begin
    I_y1 = Integral(y1 -> begin
    I_z1 = Integral(z1 -> begin
    I_x2 = Integral(x2 -> begin
    I_y2 = Integral(y2 -> begin
    I_z2 = Integral(z2 -> begin
        dx = x2 - x1; dy = y2 - y1; dz = z2 - z1;
        r2 = dx*dx + dy*dy + dz*dz
        return dx / (r2 * sqrt(r2)) end, I01)
        return integrate(I_z2, kv) end, I01)
        return integrate(I_y2, kv) end, I23)
        return integrate(I_x2, kv) end, I01)
        return integrate(I_z1, kv) end, I01)
        return integrate(I_y1, kv) end, I01)
    return integrate(I_x1, kv)
end

end # dn2