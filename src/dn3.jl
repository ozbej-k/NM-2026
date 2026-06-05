module dn3

using LinearAlgebra

export InitialProblem, ODESolution, PendulumParams
export DOPRI5, solve, pendulum_1st_order, period_from_solution, math_pendulum_energy

"""
Definition struct for describing an ODE initial problem.
"""
struct InitialProblem{TU,TT,TP}
    f # right-hand side of ODE u' = f(t, u, p)
    u0::TU # initial value
    tint::Tuple{TT,TT} # interval on which we seek the solution
    p::TP # system parameters
end

"""
Data struct for ODE solutions containing the initial problem reference, time values, their solution approximations and derivative values.
"""
struct ODESolution{TU,TT,TP}
    zp::InitialProblem{TU,TT,TP} # reference to the initial problem
    t::Vector{TT} # time values (independent variable)
    u::Vector{TU} # approximations of the solution values
    du::Vector{TU} # computed derivative values
end

"""
Parameters for the DOPRI5 method.
"""
struct DOPRI5{T}
    h::T # initial step length
    atol::T # absolute tolerance
    rtol::T # relative tolerance
end

"""
Solves initial problem for ODE using the DOPRI5 method with adaptive step control.
"""
function solve(zp::InitialProblem{TU,TT,TP}, method::DOPRI5) where {TU,TT,TP}
    t0, tk = zp.tint
    u0 = zp.u0
    h = method.h
    direction = sign(tk - t0)
    t = [t0]
    u = [u0]
    du = TU[]
    while (tk - t0) * direction > 0
        if abs(h) > abs(tk - t0)
            h = tk - t0
        end

        t_candidate, u_candidate, du0, h_new = step(method, zp.f, t0, u0, zp.p, h)
        if t_candidate != t0
            t0 = t_candidate
            u0 = u_candidate
            push!(t, t0)
            push!(u, u0)
            push!(du, du0)
        end
        h = h_new
    end

    push!(du, zp.f(t[end], u[end], zp.p))
    ODESolution(zp, t, u, du)
end


"""
Calculate a step of the DOPRI5 method.
"""
function step(m, f, t0, u0, p, h)
    # table
    k1 = f(t0, u0, p)
    k2 = f(t0 + h * (1 / 5), u0 + h * (1 / 5) * k1, p)
    k3 = f(t0 + h * (3 / 10), u0 + h * ((3 / 40) * k1 + (9 / 40) * k2), p)
    k4 = f(t0 + h * (4 / 5), u0 + h * ((44 / 45) * k1 - (56 / 15) * k2 + (32 / 9) * k3), p)
    k5 = f(t0 + h * (8 / 9), u0 + h * ((19372 / 6561) * k1 - (25360 / 2187) * k2 + (64448 / 6561) * k3 - (212 / 729) * k4), p)
    k6 = f(t0 + h, u0 + h * ((9017 / 3168) * k1 - (355 / 33) * k2 + (46732 / 5247) * k3 + (49 / 176) * k4 - (5103 / 18656) * k5), p)

    # 5th order
    u5 = u0 + h * ((35 / 384) * k1 + (500 / 1113) * k3 + (125 / 192) * k4 - (2187 / 6784) * k5 + (11 / 84) * k6)

    # 4th order
    k7 = f(t0 + h, u5, p)
    u4 = u0 + h * ((5179 / 57600) * k1 + (7571 / 16695) * k3 + (393 / 640) * k4 - (92097 / 339200) * k5 + (187 / 2100) * k6 + (1 / 40) * k7)

    error_estimate = norm(u5 - u4)
    tol = m.atol + m.rtol * max(norm(u0), norm(u5))

    # step control
    if error_estimate == 0
        h_new = 2h
        accept = true
    else
        err = error_estimate / tol
        fac = 0.9 * err^(-1 / 5)
        fac = clamp(fac, 0.2, 5.0)
        h_new = h * fac
        accept = error_estimate <= tol
    end

    if accept
        return t0 + h, u5, k7, h_new
    else # reject step
        return t0, u0, k7, h_new
    end
end

"""
Parameters for the mathematical pendulum system.
"""
struct PendulumParams
    g::Float64
    l::Float64
end

"""
Returns the first-order system for the mathematical pendulum equation.

Where u is the state vector: [angle, velocity]
"""
function pendulum_1st_order(_, u, p)
    return [u[2], -(p.g / p.l) * sin(u[1])]
end

"""
Energy calculation for the mathematical pendulum based on initial angle.
"""
math_pendulum_energy(theta_0, g, l) = g * l * (1 - cos(theta_0))

"""
Detects oscilation period from given ODESolution.
"""
function period_from_solution(sol::ODESolution)
    t = sol.t
    omega = [ui[2] for ui in sol.u]
    turning = Float64[]
    for i in 2:lastindex(omega)
        if omega[i-1] > 0 && omega[i] <= 0
            alpha = omega[i-1] / (omega[i-1] - omega[i])
            tc = t[i-1] + alpha * (t[i] - t[i-1])
            push!(turning, tc)
        end
    end
    return turning[2] - turning[1]
end

end # dn3
