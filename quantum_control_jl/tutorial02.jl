using QuantumPropagators: hamiltonian, propagate
using OrdinaryDiffEq
using Optimization, OptimizationNLopt

using Plots

# Set up thicker default lines in plots
Plots.default(
    linewidth               = 2.0,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8)
)

# Some utilities for showing hints and solutions
include("tutorial_utils/exercise_2_TLS.jl")

"""
Blackman shaped pulse

# Keyword Arguments

* `ω`: carrier frequency of the pulse
* `E₀`: pulse amplitude
* `ti`: start of the pulse
* `tf`: end of the pulse
* `α`: Blackman parameter, should be `0.16`
"""
function blackman_pulse(t; ω, E₀, ti, tf, a=0.16)
    ΔT = tf - ti
    if (t < ti) || (t > tf)
        return 0.0
    else
        return (E₀/2) * cos(ω*t) * (1 - a - cos(2π * (t-ti)/ΔT) + a * cos(4π * (t-ti)/ΔT))
    end
end


"""
Two-level-system Hamiltonian and canonical eigenstates.

# Arguments

* `ω`: Energy separation of the qubit levels
* `ϵ`: Control function
"""
function ham_and_states(ω, ϵ)
    σz = Float64[1 0; 0 -1]
    σx = Float64[0 1; 1 0]
    H₀ = -(1/2) * ω * σz
    ket0 = ComplexF64[1, 0]
    ket1 = ComplexF64[0, 1]
    return hamiltonian(H₀, (σx, ϵ)), ket0, ket1
end

T = 10.0
nt = 5001

tlist = collect(range(0, T; length=nt));

ω = 12.0  # carrier frequency, should be sufficiently large
E₀ = 0.5  # pulse amplitude
ΔT = 5.0  # pulse length
# make pulse symmetric around the middle of the time interval
ti = T/2 - ΔT/2
tf = T/2 + ΔT/2

ϵ(t) = blackman_pulse(t; ω, E₀, ti, tf)
H, ket0, ket1 = ham_and_states(ω, ϵ);

states = propagate(ket0, H, tlist; method=OrdinaryDiffEq, storage=true)

using QuantumPropagators.Controls: discretize

function plot_population_and_pulse(tlist, states, ϵ)
    E = discretize(ϵ, tlist)
    plot(tlist, abs.(E) / maximum(E); label="|E|", color=:lightgray)
    plot!(tlist, abs2.(states)'; label=["|0⟩" "|1⟩"], ls=[:solid :dash], color=["#1f77b4" "#ff7f0e"])
    plot!(; xlabel="time", ylabel="population", legend=:right)
    infidelity = 1 - abs2.(states[:,end])[2]
    annotate!(
        0, 0.9,
        ("infidelity: $(round(infidelity; digits=5))", 10, :left)
    )
end

plot_population_and_pulse(tlist, states, ϵ)


function evolve_and_plot_parameterized_pulse(; E₀, ΔT, ω=12.0, T=10, nt=5001)
    ti = T/2 - ΔT/2
    tf = T/2 + ΔT/2
    ϵ(t) = blackman_pulse(t; ω, E₀, ti, tf)
    H, ket0, ket1 = ham_and_states(ω, ϵ)
    tlist = collect(range(0, T; length=nt));
    states = propagate(ket0, H, tlist; method=OrdinaryDiffEq, storage=true)
    plot_population_and_pulse(tlist, states, ϵ)
end

# from hint
evolve_and_plot_parameterized_pulse( E₀=0.93, ΔT=8.0 )

# Try to optimize
using UnPack: @unpack  # unpack NamedTuple into variables
function loss(x, constants)
    E₀, ΔT = x
    @unpack ω, T, nt = constants
    ti = T/2 - ΔT/2
    tf = T/2 + ΔT/2
    tlist = collect(range(0, T; length=nt));
    ϵ(t) = blackman_pulse(t; ω, E₀, ti, tf)
    tlist = collect(range(0, T; length=nt));
    H, ket0, ket1 = ham_and_states(ω, ϵ)
    Ψout = propagate(ket0, H, tlist; method=OrdinaryDiffEq)
    fidelity = abs2(Ψout[2])
    return 1.0 - fidelity
end

guess = [0.5, 5.0]
upper_bounds = [10.0, T]
prob = OptimizationProblem(
    loss,
    guess,
    (; ω, T, nt),  # this is a NamedTuple, forwarding the global variables
    lb=[0.0, 0.0],
    ub=upper_bounds,
    stopval=(1-0.999),  # below which error to stop the optimization
);

loss(guess, (; ω, T, nt))

obtained_fidelities = Float64[];  # for keeping track of the fidelity in each iteration
function callback(state, loss_val)
    global obtained_fidelities
    fid = 1 - loss_val
    push!(obtained_fidelities, fid)
    print("Iteration: $(length(obtained_fidelities)), current fidelity $(round(fid; digits=4))\r")
    return false
end

obtained_fidelities = Float64[];
res = Optimization.solve(prob, NLopt.LN_NELDERMEAD(); maxiters=500, callback)
println("\n\nHighest fidelity reached: $(round((1 - res.objective) * 100; digits=2))%")
if res.objective < 1e-3
    println("\tcongratulations, you have obtained population inversion!")
else
    println("\tbad guess, please try again!")
end
plot(obtained_fidelities; marker=:cross, label="", xlabel="optimization iteration", ylabel="fidelity")

evolve_and_plot_parameterized_pulse(E₀=res.u[1], ΔT=res.u[2])
println("E₀ = $(round(res.u[1]; digits=3))\nΔT = $(round(res.u[2]; digits=3))")