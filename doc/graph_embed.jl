using dn1, Graphs, GraphRecipes, Plots

function lestev(omega)
    G = krožna_lestev(8)
    t = range(0, 2pi, 9)[1:end-1]
    x = cos.(t)
    y = sin.(t)
    točke = hcat(hcat(x, y)', zeros(2, 8))
    # funkcija hcat zloži vektorje po stolpcih v matriko
    fix = 1:8
    it = vloži!(G, fix, točke, omega)
    return it, G, točke
end

function grid(omega)
    m, n = 6, 6
    G = Graphs.grid((m, n), periodic=false)
    rob = filter(v -> degree(G, v) <= 3, vertices(G))
    urejen_rob = [rob[1]]
    # uredi točke na robu v cikel
    for i = 1:length(rob)-1
        sosedi = neighbors(G, urejen_rob[end])
        sosedi = intersect(sosedi, rob)
        sosedi = setdiff(sosedi, urejen_rob)
        push!(urejen_rob, sosedi[1])
    end
    t = range(0, 2pi, length(rob) + 1)[1:end-1]
    točke = zeros(2, n * m)
    točke[:, urejen_rob] = hcat(cos.(t), sin.(t))'
    it = vloži!(G, urejen_rob, točke, omega)
    return it, G, točke
end

omegas = 0.1:0.01:1.9
iterations = [lestev(omega)[1] for omega in omegas]

println("Plotting unembedded cyclic ladder")
G = krožna_lestev(8)
display(graphplot(G, curves=false))
println("Press enter to continue")
readline()

println("Plotting embedded cyclic ladder")
best_omega = omegas[argmin(iterations)]
_, G, točke = lestev(best_omega)
display(graphplot(G, x=točke[1, :], y=točke[2, :], curves=false))
println("Press enter to continue")
readline()

println("Best omega for ladder: $best_omega")
println("Plotting SOR iterations needed per omega choice")

display(plot(omegas, iterations, xlabel="\\omega", ylabel="SOR iterations"; legend=false))
println("Press enter to continue")
readline()


omegas = 0.1:0.01:1.9
iterations = [grid(omega)[1] for omega in omegas]

println("Plotting unembedded grid")
G = Graphs.grid((6, 6), periodic=false)
display(graphplot(G, curves=false))
println("Press enter to continue")
readline()

println("Plotting embedded grid")
best_omega = omegas[argmin(iterations)]
_, G, točke = grid(best_omega)
display(graphplot(G, x=točke[1, :], y=točke[2, :], curves=false))
println("Press enter to continue")
readline()

println("Plotting SOR iterations needed per omega choice")
println("Best omega for grid: $best_omega")

display(plot(omegas, iterations, xlabel="\\omega", ylabel="of SOR iterations"; legend=false))

