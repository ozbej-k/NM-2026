using dn3, Plots

# plot mathematical vs harmonic pendulum oscilation
println("Plotting mathematical vs harmonic pendulum oscilation")
g = 9.81; l = 1.0
theta0 = pi/4
p = PendulumParams(g,l)
problem = InitialProblem(pendulum_1st_order, [theta0,0.0], (0.0,20.0), p)
sol = solve(problem, DOPRI5(1e-3, 1e-13, 1e-13))

function harmonic_theta(t, theta0, thetadot0, g, l)
    omega0 = sqrt(g / l)
    return theta0 * cos(omega0 * t) + (thetadot0 / omega0) * sin(omega0 * t)
end

theta_math = [ui[1] for ui in sol.u]
theta_harm = [harmonic_theta(t, theta0, 0.0, g, l) for t in sol.t]

plt = plot(sol.t, theta_math, lw = 2, label = "Mathematical")
plot!(sol.t, theta_harm, lw = 2, linestyle = :dash, label = "Harmonic")
display(plt)
println("Press enter to continue")
readline()

# plot mathematical vs harmonic pendulum energy
println("Plotting mathematical vs harmonic pendulum energy")
energies = Float64[]; periods = Float64[]
for theta0 in range(0.01, pi - 0.1, length=20)
    _p = PendulumParams(g, l)
    _problem = InitialProblem(pendulum_1st_order, [theta0, 0.0], (0.0, 20.0), _p)
    _sol = solve(_problem, DOPRI5(1e-3, 1e-13, 1e-13))
    push!(energies, math_pendulum_energy(theta0, g, l))
    push!(periods, period_from_solution(_sol))
end

harmonic_period(g, l) = 2pi * sqrt(l / g)

plot(energies, periods, xlabel="Energy", ylabel="Oscilation period", lw=2, grid=true)
hline!([harmonic_period(g, l)], linestyle=:dash, label="Harmonic period", lw=2)    


