using dn3, DifferentialEquations, LinearAlgebra, Test

@testset "DOPRI5 accuracy" begin
    g = 9.81; l = 1.0
    p = PendulumParams(g, l)
    θ0 = 1.0; ω0 = 0.0
    tspan = (0.0, 10.0)

    # my solver
    prob_my = InitialProblem(pendulum_1st_order, [θ0, ω0], tspan, p)
    sol_my = dn3.solve(prob_my, DOPRI5(1e-3, 1e-13, 1e-13))

    # Julia solver
    function f!(du, u, p, t)
        du[1] = u[2]
        du[2] = -(p.g / p.l) * sin(u[1])
    end
    prob_ref = ODEProblem(f!, [θ0, ω0], tspan, p)
    sol_ref = DifferentialEquations.solve(prob_ref, Tsit5(), abstol=1e-13, reltol=1e-13)

    # compare solvers
    for (t, u_my) in zip(sol_my.t, sol_my.u)
        u_ref = sol_ref(t)
        @test abs(u_my[1] - u_ref[1]) < 1e-10
        @test abs(u_my[2] - u_ref[2]) < 1e-10
    end
end

