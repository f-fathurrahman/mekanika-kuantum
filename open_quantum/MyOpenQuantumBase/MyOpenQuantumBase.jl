module MyOpenQuantumBase

using DocStringExtensions

import LinearAlgebra: mul!, lmul!, axpy!, ishermitian, Hermitian, eigen, eigen!,
    tr, diag, Diagonal, norm, I, Bidiagonal, qr, diagm
import StaticArrays: SMatrix, MMatrix
import SparseArrays: sparse, issparse, spzeros, SparseMatrixCSC
import LinearAlgebra.BLAS: her!, gemm!
import QuadGK: quadgk!, quadgk

"""
$(TYPEDEF)

This is a supertype for Hamiltonians with elements of type `T`. Any object of a
subtype of `AbstractHamiltonian` should implement two methods corresponding to 
the call signatures `H(t)` and `H(du, u, p, t)`.

`H(t)`: This method, when called with a single argument `t`, should return the 
value of the Hamiltonian at time `t`.

`H(du, u, p, t)`: This method should update the `du` object in-place according to 
the evolution equation ``du += -i H(t) u`` when ``u`` is a vector, or 
``du += -i [H(t), u]`` when ``u`` is a matrix. Here, `u` is the state vector or 
density matrix, `p` contains parameters for the ODE solver, and `t` is the 
current time.
"""
abstract type AbstractHamiltonian{T<:Number} end

"""
$(TYPEDEF)

Base for types defining Hamiltonians using dense matrices.
"""
abstract type AbstractDenseHamiltonian{T<:Number} <: AbstractHamiltonian{T} end


"""
$(TYPEDEF)

Base for types defining Hamiltonians using sparse matrices.
"""
abstract type AbstractSparseHamiltonian{T<:Number} <: AbstractHamiltonian{T} end


"""
$(TYPEDEF)

Base for types defining quantum annealing process.
"""
abstract type AbstractAnnealing{constant_hamiltonian} end

"""
$(TYPEDEF)

Base for types defining system bath coupling operators in open quantum system models.
"""
abstract type AbstractCouplings end

"""
$(TYPEDEF)

Base for types defining bath object.
"""
abstract type AbstractBath end

"""
$(TYPEDEF)
"""
abstract type AbstractLiouvillian end

include("lobpcg.jl")
include("base_util.jl")
include("unit_util.jl")
include("math_util.jl")
include("matrix_util.jl")
include("interpolation.jl")
include("dev_tools.jl")

include("coupling/coupling.jl")
include("coupling/interaction.jl")

include("bath/bath_util.jl")
include("bath/ohmic.jl")
include("bath/hybrid_ohmic.jl")
include("bath/spin_fluc.jl")
include("bath/custom.jl")

include("integration/cpvagk.jl")
include("integration/ext_util.jl")

include("hamiltonian/hamiltonian_base.jl")
include("hamiltonian/interface.jl")
include("hamiltonian/dense_hamiltonian.jl")
include("hamiltonian/sparse_hamiltonian.jl")
include("hamiltonian/static_array_hamiltonian.jl")
include("hamiltonian/adiabatic_frame_hamiltonian.jl")
include("hamiltonian/interp_hamiltonian.jl")
include("hamiltonian/custom_hamiltonian.jl")

include("annealing/annealing_type.jl")
include("annealing/displays.jl")

include("opensys/opensys_util.jl")
include("opensys/diffeq_liouvillian.jl")
include("opensys/redfield.jl")
include("opensys/davies.jl")
include("opensys/cgme.jl")
include("opensys/lindblad.jl")
include("opensys/stochastic.jl")
include("opensys/trajectory_jump.jl")

include("projection/projection.jl")

export SparseHamiltonian, DenseHamiltonian, AdiabaticFrameHamiltonian,
    InterpDenseHamiltonian, InterpSparseHamiltonian, CustomDenseHamiltonian, Hamiltonian, ConstantHamiltonian
export rotate, isconstant, isdimensionlesstime
export eigen_decomp, haml_eigs

export temperature_2_β, temperature_2_freq, β_2_temperature, freq_2_temperature
export σx, σz, σy, σi, σ, ⊗, PauliVec, spσx, spσz, spσi, spσy, σ₊, σ₋,
    bloch_to_state, creation_operator, annihilation_operator, number_operator

export q_translate, standard_driver, local_field_term, two_local_term,
    single_clause, q_translate_state, collective_operator, hamming_weight_operator

export check_positivity, check_unitary, check_density_matrix, partial_trace,
    matrix_decompose, low_level_matrix, fidelity, inst_population, gibbs_state, purity, check_pure_state, find_degenerate

export haar_unitary

export construct_interpolations, gradient, log_uniform

export hamiltonian_from_function, evaluate, issparse, get_cache, update_cache!,
    update_vectorized_cache!

export ConstantCouplings, TimeDependentCoupling, TimeDependentCouplings,
    CustomCouplings, collective_coupling, Interaction, InteractionSet

export Annealing, ODEParams, set_u0!, Evolution

export Lindblad, EnsembleFluctuator, DiffEqLiouvillian

export Ohmic, OhmicBath, CustomBath, CorrelatedBath, ULEBath, HybridOhmic,
    HybridOhmicBath
export correlation, polaron_correlation, γ, S, spectrum, info_freq

export τ_B, τ_SB, coarse_grain_timescale

export InplaceUnitary, EᵨEnsemble, sample_state_vector

export ProjectedSystem, project_to_lowlevel, get_dθ, concatenate

# APIs for test and developement tools, may move to a different repo latter
export random_ising, alt_sec_chain, build_example_hamiltonian

end  # module MyOpenQuantumBase
