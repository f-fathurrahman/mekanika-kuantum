Pkg.activate("quantum_goatjl", shared=true)
push!(LOAD_PATH, "./GOAT")

using GOAT
using LinearAlgebra, SparseArrays
using OrdinaryDiffEq
using Optim, LineSearches

σx = sparse(ComplexF64[0 1 ; 1  0])# The Pauli X matrix
σy = sparse(ComplexF64[0 -im ; im 0]) # The Pauli Y matrix
σz = sparse(ComplexF64[1 0 ; 0 -1]) # The Pauli Z matrix
σI = sparse(ComplexF64[1 0 ; 0  1]) # The identity matrix

# Hamiltonian parameters and the control basis operators
ω = 5*2pi # Setting a qubit frequency
T = 10.0 # in arbitrary time units
H0 = (ω/2)*σz # Setting the drift component of the Hamiltonian
basis_ops = [σx] # The control basis operator

# Control ansatz, one basis only
N_basis_funcs = 1
# The control ansatz is a gaussian envelope of a cosine wave at frequency ω
Ω(p,t,i) = gaussian_coefficient(p,t,i,N_basis_funcs)*cos(ω*t)
# amplitude, center, width

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
ODE_options = (abstol=1e-9, reltol=1e-9, alg=Vern9())

# Define the optimizer options from Optim.jl (See Optim.jl docs for info)
optim_alg = Optim.LBFGS(linesearch=LineSearches.BackTracking()) # A Back-Tracking linesearch
optim_options = Optim.Options(
    g_tol=1e-12,
    iterations=100,
    store_trace=true,
    show_trace=true,
    extended_trace=false,
    allow_f_increases=false
)

p0 = [0.5, T/2, T/8] # The initial parameter guesses
opt_param_inds = [1]
# The parameters of the vector p0 to optimize
# (just the amplitude parameter -- p0[1] -- in this case)

derivs_per_core = 1

# Next we put all of these parameters into a `QOCParameters` struct:

params = QOCParameters(
    ODE_options,
    SE_reduce_map,
    GOAT_reduce_map,
    optim_alg,
    optim_options;
    derivs_per_core=derivs_per_core
)

# Finally we run our optimization
res = find_optimal_controls(p0, opt_param_inds, sys, prob, params)

using Plots, PlotThemes
theme(:dark)

t_grid = collect(0:0.01:T)
control_signal = similar(t_grid)
p_opt = copy(p0)
for (idx,v) in enumerate(opt_param_inds)
    p_opt[v] = Optim.minimizer(res)[idx]
end
Npoints = length(t_grid)
for i in 1:Npoints
    control_signal[i] = Ω(p_opt, t_grid[i], 1)
end
fig_control = plot(t_grid, control_signal)

# Solve
sol = solve_SE(sys, T, p_opt) # The solutions are unitary gates
ψ0 = complex.([1.0, 0.0])
ψ = copy(ψ0)
U = zeros(ComplexF64, 2, 2)
pop0 = zeros(Float64, Npoints)
pop1 = zeros(Float64, Npoints)
# Propagate using
for i in 1:Npoints
    t = t_grid[i]
    U[:,:] = sol(t) # interpolate
    ψ[:] = U * ψ
    c2 = real(ψ' * ψ)
    ψ[:] = ψ/sqrt(c2)
    pop0[i] = abs2(ψ[1])
    pop1[i] = abs2(ψ[2])
    #println("t = $t ψ = $(ψ)")
end
fig_pop = plot(t_grid, pop0)
plot!(fig_pop, t_grid, pop1)
