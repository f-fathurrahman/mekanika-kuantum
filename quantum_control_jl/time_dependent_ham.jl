struct TimeDepHamiltonian
    H0::Matrix{ComplexF64}
    Hs::Union{Vector{Matrix{ComplexF64}},Nothing}
    amplitudes::Union{Vector{Function},Nothing}
end

function evaluate_ham(Ham::TimeDepHamiltonian, t::Float64)
    Ht = similar(Ham.H0)
    evaluate_ham!(Ham, t, Ht)
    return Ht
end

function evaluate_ham!(
    Ham::TimeDepHamiltonian,
    t::Float64,
    Ht::Matrix{ComplexF64}
)
    H0 = Ham.H0
    Hs = Ham.Hs
    amplitudes = Ham.amplitudes
    #
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



# Time-discretized Hamiltonian (uniform time grid)
# The amplitudes are discretized, using the same discret time grid
struct DiscreteTimeDepHamiltonian
    H0::Matrix{ComplexF64}
    Hs::Union{Vector{Matrix{ComplexF64}},Nothing}
    amplitudes::Union{Vector{Vector{Float64}},Nothing}
end


function evaluate_ham!(
    Ham::DiscreteTimeDepHamiltonian,
    idx_t::Int64,
    Ht::Matrix{ComplexF64}
)
    H0 = Ham.H0
    Hs = Ham.Hs
    amplitudes = Ham.amplitudes
    #
    @views Ht[:] = H0[:]
    # Early return
    if isnothing(Hs)
        return
    end
    Nterms = size(Hs, 1)
    for i in 1:Nterms
        Ht[:,:] .+= amplitudes[i][idx_t] .* Hs[i][:,:]
    end
    return
end
