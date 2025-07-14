using Infiltrator

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

function control_function(t)
    t0 = 25.0
    ωl = 10.0
    E0 = 1.2
    ϕ = 0
    τ = 1.5
    return E0 * cos(ωl * (t - t0) + ϕ) * exp(-(t - t0)^2 / (2τ^2));
end

function main_time_dep_01()


    ω = 10.0
    μ01 = 1.0
    H0 = -0.5 * ω * ComplexF64[
        1   0
        0  -1
    ]
    H1 = -μ01 * ComplexF64[
        0  1
        1  0
    ]

    Ham = TimeDependentHamiltonian(H0, [H1], [control_function])

    tlist = collect(range(0, 50.0, length=10_000))
    pulse = control_function.(tlist)
    
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

    @exfiltrate

end