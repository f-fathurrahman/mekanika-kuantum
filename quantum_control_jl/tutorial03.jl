using QuantumControl
using QuantumPropagators: propagate
using QuantumPropagators: ExpProp

using Krotov

using Plots

# Set up thicker default lines in plots
Plots.default(
    linewidth               = 2.0,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8)
)

"""Two-level-system Hamiltonian

# Keyword Arguments

* `omega` (float): energy separation of the qubit levels
* `eps0` (function): amplitude eps0(t) of the driving field
"""
function ham_and_states(; omega=1.0, eps0=(t -> 1.0))
    H₀ = -0.5 * omega * [
        1   0
        0  -1
    ]
    H₁ = Float64[
        0  1
        1  0
    ]
    Ψ₀ = ComplexF64[1, 0] # State |0⟩
    Ψ₁ = ComplexF64[0, 1] # State |1⟩
    H = hamiltonian(H₀, (H₁, eps0))
    return H, Ψ₀, Ψ₁
end

"""Shape function for the field update"""
S(t) = QuantumControl.Shapes.flattop(t; T=10.0, t_rise=0.5, func=:sinsq);

tlist = collect(range(0, 10; length=81));

"""Guess Amplitude (unshaped)"""
E(t; A=0.1, σ=2) = A * exp(-(t-5)^2 / (2 * σ^2)) * cos(3t);

H, Ψ₀, Ψ₁ = ham_and_states(eps0=(t -> S(t) * E(t)));

using QuantumControl.Controls: discretize

function plot_pulse(pulse, tlist)
    fig = plot(; xlabel="time", ylabel="pulse amplitude")
    plot!(fig, tlist, discretize(pulse, tlist); label="")
    return fig
end
plot_pulse(H.amplitudes[1], tlist)

states = propagate(Ψ₀, H, tlist; method=ExpProp, storage=true);
plot(abs2.(states)', labels=["0" "1"]; xlabel="time", ylabel="population")
plot(abs2.(states[2,:]), label="1"; xlabel="time", ylabel="population")

const 𝕚 = 1im
σ_x = ComplexF64[0 1; 1 0]
σ_y = ComplexF64[0 -𝕚; 𝕚 0]
σ_z = ComplexF64[1 0; 0 -1]

bloch_vals = propagate(Ψ₀, H, tlist; method=ExpProp, observables=[σ_x, σ_y, σ_z], storage=true);

using QuantumControl.Functionals: J_T_ss
trajectories = [Trajectory(Ψ₀, H; target_state=Ψ₁)];

problem = ControlProblem(
    trajectories,
    tlist;
    prop_method=ExpProp,
    J_T=J_T_ss,
    iter_stop=50,
);

res = optimize(problem; method=Krotov, lambda_a=25, update_shape=S)

using QuantumPropagators.Controls: get_controls, substitute
H_opt = substitute(H, Dict(get_controls(H)[1] => res.optimized_controls[1]));
states_opt = propagate(Ψ₀, H_opt, tlist; method=ExpProp, storage=true);
plot_pulse(H_opt.amplitudes[1], tlist)
plot(abs2.(states_opt)', labels=["0" "1"]; xlabel="time", ylabel="population", legend=:right)
bloch_vals = propagate(Ψ₀, H_opt, tlist; method=ExpProp, observables=[σ_x, σ_y, σ_z], storage=true);

#plot_bloch(bloch_vals)
Ntrajs = size(bloch_vals, 2)
bloch_trajs = Vector{Vector{Float64}}(undef,Ntrajs)
for i in 1:Ntrajs
    bloch_trajs[i] = bloch_vals[:,i]
end
