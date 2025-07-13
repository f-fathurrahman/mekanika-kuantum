using Test
using MyQuantumPropagators
using MyQuantumPropagators: Cheby
using MyQuantumPropagators.Storage
using LinearAlgebra

# We use Rabi cycling in a TLS as a test case, which allows to to compare
# the propagation with the known analytic solution.

function debug_main()

    Ψ0 = ComplexF64[1, 0]
    Ĥ = ComplexF64[
         0   0.5
        0.5   0
    ]
    tlist = collect(range(0, 1.5π, length=101)) # 3π/2 pulse

    generator = (Ĥ,)

    storage = init_storage(Ψ0, tlist)
    @test isa(storage, Matrix)
    @test eltype(storage) == ComplexF64

    Ψ_out = propagate(Ψ0, generator, tlist; method=:expprop, storage=storage)
    Ψ_expected = ComplexF64[-1/√2, -1im/√2]  # note the phases

    pop0 = abs.(storage[1, :]) .^ 2
    #println(lineplot(tlist ./ π, pop0, ylim=[0, 1], title="fw prop"))

    @test norm(Ψ_out - Ψ_expected) < 1e-12
    @test pop0[end] ≈ 0.5

#=

    # Propagating backward in time should exactly reverse the dynamics (since
    # they are unitary). Thus, we should end up back at the initial state, and
    # the stored states (being filled back-to-front) should exactly match the
    # stored states from the forward propagation.

    storage_bw = init_storage(Ψ0, tlist)
    Ψ_out_bw = propagate(
        Ψ_out,
        generator,
        tlist;
        method=:expprop,
        backward=true,
        storage=storage_bw
    )

    pop0_bw = abs.(storage_bw[1, :]) .^ 2
    #println(lineplot(tlist ./ π, pop0_bw, ylim=[0, 1], title="bw prop"))

    @test norm(Ψ_out_bw - Ψ0) < 1e-12
    @test pop0_bw[1] ≈ 1.0
    @test norm(storage - storage_bw) < 1e-12
=#


end

