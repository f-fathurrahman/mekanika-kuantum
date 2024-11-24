# from QuantumDynamics.jl paper

using QuantumDynamics

H = [-0.1im -1.0; -1.0 -0.5im]
V(t) = 12*cos(10*t)
EF = Utilities.ExternalField(V, [1.0+0.0im 0.0; 0.0 -1.0])

# define the initial condition, the time step and the
# number of steps of simulation
ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
dt = 0.125
ntimes = 100

times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes, external_fields=[EF])
