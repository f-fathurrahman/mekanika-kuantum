using Infiltrator 
using MyQuantumPropagators
using LinearAlgebra

function debug_main()
    psi0 = ComplexF64[1, 0]
    H0 = ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0, 0.1π, length=2))
    generator = (H0,)
    psi1 = propagate(psi0, generator, tlist; method=:expprop, storage=true)
    @exfiltrate

    #= Compare with
    Δt = tlist[2] - tlist[1]
    psi_new = exp(-1m * H0 * Δt) * psi0
    =#
end

