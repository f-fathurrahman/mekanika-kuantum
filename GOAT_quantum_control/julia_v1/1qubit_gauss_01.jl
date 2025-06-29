Pkg.activate("quantum_goatjl", shared=true)

using Revise
using Infiltrator

using Distributed: pmap
using LinearAlgebra
using SparseArrays

#using LinearAlgebra, SparseArrays # For instantiating necessary matrices
using OrdinaryDiffEq # Load differetial equation algorithms
using Optim, LineSearches # Load optimization and linesearch algorithms

includet("goat_functions.jl")

function main_debug()

    ## Defining the `ControllableSystem`

    # We begin by defining the necessary operators for the Hamiltonian
    σx = sparse(ComplexF64[0 1 ; 1  0])# The Pauli X matrix
    σy = sparse(ComplexF64[0 -im ; im 0]) # The Pauli Y matrix
    σz = sparse(ComplexF64[1 0 ; 0 -1]) # The Pauli Z matrix
    σI = sparse(ComplexF64[1 0 ; 0  1]) # The identity matrix

    # Next we define the Hamiltonian parameters and the control basis operators

    ω = 5*2pi # Setting a qubit frequency
    T = 10.0 # in arbitrary time units
    H0 = (ω/2)*σz # Setting the drift component of the Hamiltonian
    basis_ops = [σx] # The control basis operator

    # Next we choose a control ansatz specifying only a single basis function

    N_basis_funcs = 1
    # The control ansatz is a gaussian envelope of a cosine wave at frequency ω
    Ω(p,t,i) = gaussian_coefficient(p,t,i,N_basis_funcs)*cos(ω*t)

    # The derivative of the control ansatz
    ∂Ω(p,t,i,l) = ∂gaussian_coefficient(p,t,i,l,N_basis_funcs)*cos(ω*t)

    # Moving into the interaction picture
    rotating_frame_generator = -Diagonal(Matrix(H0))

    sys = ControllableSystem(H0, basis_ops, rotating_frame_generator, Ω, ∂Ω)
    
    # First we define our target unitary operator
    U_target = copy(σx) # The targt unitary: X|0> = |1>

    # Now, we provide a projector onto the computational and ancillary subspaces

    # The projector onto the computational subspace (which is the whole Hilbert space in this case)
    Pc = σI
    # The projector onto the ancillary subspaces (which is a null operator in this case)
    Pa = 0.0*Pc

    # The quantum optimal control problem. 
    prob = QOCProblem(U_target, T, Pc, Pa)

    # reduction map
    SE_reduce_map = SE_infidelity_reduce_map 
    GOAT_reduce_map = GOAT_infidelity_reduce_map


    # Define options for DifferentialEquations.jl (see DifferentialEquations.jl docs for info)
    ODE_options = (abstol = 1e-9, reltol= 1e-9, alg=Vern9())

    # Define the optimizer options from Optim.jl (See Optim.jl docs for info)
    optim_alg = Optim.LBFGS(linesearch=LineSearches.BackTracking()) # A Back-Tracking linesearch
    optim_options = Optim.Options(g_tol=1e-12,
                            iterations=10,
                            store_trace=true,
                            show_trace=true, extended_trace=false, allow_f_increases=false)

    p0 = [0.5,T/2,T/8] # The initial parameter guesses
    opt_param_inds = [1]
    # The parameters of the vector p0 to optimize (just the amplitude parameter -- p0[1] -- in this case)
    
    derivs_per_core = 1
    # This `derivs_per_core` variable specifies how many derivatives are propogated in
    # each GOAT EOMs and informs parallelization. 

    # For example, if `derivs_per_core=5`and there are 5 total parameters, then no
    # parallelization is performed. In contrast, if `derivs_per_core=2` and there are
    # 5 total parameters, then 3 processes are run in parallel: the first processor computes
    # the EOMs for 2 parameters, the second process computes the EOMs for 2 parameters, and
    # the third computes the EOMs for 1 parameter. 

    # Next we put all of these parameters into a `QOCParameters` struct:

    params = QOCParameters(ODE_options,SE_reduce_map,GOAT_reduce_map,optim_alg,optim_options; derivs_per_core=derivs_per_core)

    # Finally we run our optimization

    res = find_optimal_controls(p0, opt_param_inds, sys, prob, params)


    @infiltrate

end

