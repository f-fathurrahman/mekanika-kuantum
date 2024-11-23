using OrdinaryDiffEq

f(u, p, t) = 1.01 * p * u
u0 = 1 / 2
p = 1.0 # ODE parameter
tspan = (0.0, 1.0)

prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
sol = solve(prob, Tsit5(), adaptive=false, dt = 0.01, reltol = 1e-8, abstol = 1e-8)