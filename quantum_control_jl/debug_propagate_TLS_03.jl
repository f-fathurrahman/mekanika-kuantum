using Infiltrator 
using MyQuantumPropagators
using LinearAlgebra

function debug_propagate_03()
    psi0 = ComplexF64[1, 0]
    H0 = ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0, 1.5π, length=101)) # 3π/2 pulse
    generator = (H0, ) # only drift Hamiltonian
    psi_ref = propagate(psi0, generator, tlist; method=:expprop, storage=true)

    Ntimes = length(tlist)
    psi_t = zeros(ComplexF64, 2, Ntimes)
    psi_t[:,1] = psi0 # set initial condition
    for i in 1:Ntimes-1
        Δt = tlist[i+1] - tlist[i]
        # H0 is constant in this case
        psi_t[:,i+1] = exp(-im * H0 * Δt) * psi_t[:,i]
    end
    @exfiltrate

end

