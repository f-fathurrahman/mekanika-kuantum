using LinearAlgebra
using OrdinaryDiffEq
using QuantumPropagators

using Plots, PlotThemes
theme(:dark)
Plots.default(
    linewidth = 2.0
)

t_start = 0.0
t_stop = 50.0
Nt = 10_000
t = collect(range(t_start, t_stop, length=Nt))

ω = 10.0
μ₀₁ = 1.0
H₀ = -ω/2 * [1  0; 0 -1]

V = -μ₀₁ * [0 1; 1 0]

# Define electric field
t₀ = 25.0
ωₗ = ω/2
E₀ = 1.2
ϕ = 0
τ = 1.5
E(t) = E₀ * cos(ωₗ * (t - t₀) + ϕ) * exp(-(t - t₀)^2 / (2τ^2));

H = hamiltonian(H₀, (V, E))
Ψ₀ = ComplexF64[1, 0] # State |0⟩
Ψ₁ = ComplexF64[0, 1] # State |1⟩

# Projectors for observables
P₀ = Ψ₀ * Ψ₀'
P₁ = Ψ₁ * Ψ₁'

output = propagate(Ψ₀, H, t; method=OrdinaryDiffEq, observables=[P₀, P₁], storage=true)

E_max = maximum(abs.(E.(t)))
plot(t, abs.(E.(t)) / E_max, label="|E|")
plot!(t, real.(output[1,:]), label="|0⟩")
plot!(t, real.(output[2,:]), linestyle=:dash, label="|1⟩")
plot!(; xlabel="Time", ylabel="Population")
savefig("IMG_01.png")
