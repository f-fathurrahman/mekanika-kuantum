using Infiltrator

using MyQuantumControl
using MyQuantumPropagators: propagate
using MyQuantumPropagators: ExpProp

using MyQuantumControl.Controls: discretize
using MyQuantumControl.Functionals: J_T_ss
using MyQuantumPropagators.Controls: get_controls, substitute

using Plots, PlotThemes
theme(:dark)

struct TimeDependentHamiltonian
    H0::Matrix{ComplexF64}
    Hs::Union{Vector{Matrix{ComplexF64}},Nothing}
    amplitudes::Union{Vector{Function},Nothing}
end

function evaluate_ham(Ham::TimeDependentHamiltonian, t::Float64)
    H0 = Ham.H0
    Hs = Ham.Hs
    amplitudes = Ham.amplitudes

    N = size(H0, 1)
    # Early return
    if isnothing(Hs)
        return H0
    end
    Nterms = size(Hs, 1)
    H = copy(H0)
    for i in 1:Nterms
        H += amplitudes[i](t) * Hs[i]
    end
    return H
end

function evaluate_ham!(
    Ham::TimeDependentHamiltonian,
    t::Float64,
    Ht::Matrix{ComplexF64}
)
    H0 = Ham.H0
    Hs = Ham.Hs
    amplitudes = Ham.amplitudes

    N = size(H0, 1)
    @views Ht[:] = H0[:]
    # Early return
    if isnothing(Hs)
        return
    end
    Nterms = size(Hs, 1)
    for i in 1:Nterms
        Ht[:,:] .+= amplitudes[i](t) .* Hs[i][:,:]
    end
    return
end




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

function plot_pulse(pulse, tlist)
    fig = plot(; xlabel="time", ylabel="pulse amplitude")
    plot!(fig, tlist, discretize(pulse, tlist); label="")
    return fig
end


function debug_krotov_02()

    """Shape function for the field update"""
    S(t) = MyQuantumControl.Shapes.flattop(t; T=10.0, t_rise=0.5, func=:sinsq);

    """Guess Amplitude (unshaped)"""
    E(t; A=0.1, σ=2) = A * exp(-(t-5)^2 / (2 * σ^2)) * cos(3t);

    ω = 1.0
    H₀ = -0.5 * ω * ComplexF64[
        1   0
        0  -1
    ]
    H₁ = ComplexF64[
        0  1
        1  0
    ]

    Ham = TimeDependentHamiltonian(H₀, [H₁], [t->S(t)*E(t)])

    tlist = collect(range(0, 1.5π, length=101)) # 3π/2 pulse
    
    ket0 = ComplexF64[1.0, 0.0]
    ket1 = ComplexF64[0.0, 1.0]
    P0 = ket0 * ket0'
    P1 = ket1 * ket1'
    psi0 = ComplexF64[1.0, 0.0]
    Ntimes = length(tlist)
    psi_t = zeros(ComplexF64, 2, Ntimes)
    psi_t[:,1] = psi0 # set initial condition
    Ht = zeros(ComplexF64, 2, 2)
    for i in 1:Ntimes-1
        Δt = tlist[i+1] - tlist[i]
        evaluate_ham!(Ham, tlist[i], Ht)
        psi_t[:,i+1] = exp(-im * Ht * Δt) * psi_t[:,i]
    end

    pop0 = zeros(Float64, Ntimes)
    pop1 = zeros(Float64, Ntimes)
    for i in 1:Ntimes
        @views ψ = psi_t[:,i]
        pop0[i] = ψ' * (P0 * ψ)
        pop1[i] = ψ' * (P1 * ψ)
    end

    @infiltrate

end