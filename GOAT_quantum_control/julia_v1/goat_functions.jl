"""
    gaussian_coefficient(p,t,i,N)

Compute the time dependent coefficient for a control basis operator given by a gaussian ansatze.

# Arguments
- `p`: The vector of control parameters.
- `t`: The time at which to evaluate the coefficient function.
- `i`: The index of the control basis operator.
- `N`: The total number of basis functions used to define the coefficient function. 
"""
function gaussian_coefficient(p, t, i, N::Int64)
    c = 0.0
    K = 3
    for n = 1:N
        j = (i - 1) * K * N + (n - 1) * K # Get the linear index of the ith's control term's nth basis function
        c += p[j+1] * exp(-0.5 * ((t - p[j+2]) / p[j+3])^2)
    end
    return c
end


"""
    ∂gaussian_coefficient(p,t,i,l,N::Int64)

Computes the partial derivative of `fourier_coefficient` w.r.t p[l] evaluated at (p,t).

# Arguments
- `p`: The vector of control parameters.
- `t`: The time at which to evaluate the coefficient function.
- `i`: The index of the control basis operator.
- `l`: The index of the `p` with which to compute the partial derivative to. 
- `N`: The total number of basis functions used to define the coefficient function. 
"""
function ∂gaussian_coefficient(p, t, i, l, N::Int64)
    K = 3
    jmin = (i - 1) * K * N + 1
    jmax = i * K * N
    if l >= jmin && l <= jmax # Linear indices between the ith control term and the i+1th control terms
        c = 0.0
        for n = 1:N
            j = (i - 1) * K * N + (n - 1) * K # Get the linear index of the ith's control term's nth basis function
            if j + 1 == l
                c += exp(-0.5 * ((t - p[j+2]) / p[j+3])^2)
            elseif j + 2 == l
                c +=
                    p[j+1] * (t - p[j+2]) * exp(-0.5 * ((t - p[j+2]) / p[j+3])^2) /
                    (p[j+3]^2)
            elseif j + 3 == l
                c +=
                    p[j+1] * ((t - p[j+2])^2) * exp(-0.5 * ((t - p[j+2]) / p[j+3])^2) /
                    (p[j+3]^3)
            end
        end
        return c
    end
    return 0.0
end



"""
Instantiate a ControllableSystem struct.
"""
struct ControllableSystem{A,B,C,D,E}
    d_ms::A
    d_ls::A
    d_vs::B
    c_ls::Vector{Vector{Int64}}
    c_ms::Vector{Vector{Int64}}
    c_vs::Vector{Vector{ComplexF64}}
    coefficient_func::C
    ∂coefficient_func::D
    reference_frame_generator::E
    reference_frame_storage::E
    use_rotating_frame::Bool
    dim::Int64
    orthogonal_basis::Bool
    tol::Float64
end


"""
    ControllableSystem(drift_op, basis_ops, RF_generator::Eigen, c_func, ∂c_func; <keyword arguments>)

Instantiate a ControllableSystem struct with a specified `RF_generator`. 

The `RF_generator` is a specification of the time-independent generator of the rotating
reference frame. i.e., if ``∂V(t)=A*V(t)`` where A is the generator of type ``Eigen``.

# Keyword Arguments
- `sparse_tol::Float`: A tolerance which defines a threshold to discard matrix elements
"""
function ControllableSystem(
    drift_op,
    basis_ops,
    RF_generator::Eigen,
    c_func,
    ∂c_func;
    sparse_tol = 1e-12,
)
    d = size(drift_op, 1)
    F = RF_generator
    N = size(basis_ops, 1)
    c_ls = []
    c_ms = []
    c_vs = []
    as = F.values
    a_diffs = zeros(ComplexF64, d, d)
    aj_drift_aks = zeros(ComplexF64, d, d)
    aj_hi_aks = zeros(ComplexF64, d, N, d)
    for j = 1:d
        for k = 1:d
            aj = F.values[j]
            ak = F.values[k]
            a_diffs[j, k] = aj - ak
            aj_vec = @view F.vectors[:, j]
            ak_vec = @view F.vectors[:, k]
            aj_drift_ak = adjoint(aj_vec) * drift_op * ak_vec
            aj_drift_aks[j, k] = aj_drift_ak

            new_basis_op = sparse(aj_vec * adjoint(ak_vec))
            droptol!(new_basis_op, sparse_tol)
            for i = 1:N
                aj_hi_ak = adjoint(aj_vec) * basis_ops[i] * ak_vec
                aj_hi_aks[j, i, k] = aj_hi_ak
            end

            ls, ms, vs = findnz(new_basis_op)
            push!(c_ls, ls)
            push!(c_ms, ms)
            push!(c_vs, vs)
        end
    end

    function new_c_func(p, t, j, k)
        c = 0.0 + 0.0im
        diag_term = 0.0 + 0.0im
        if j == k
            diag_term = as[j]
        end
        for i = 1:N
            c += c_func(p, t, i) * aj_hi_aks[j, i, k]
        end
        adiff = a_diffs[j, k]
        aj_drift_ak = aj_drift_aks[j, k]
        return cis(t * adiff) * (aj_drift_ak + c + diag_term)
    end

    function new_∂c_func(p, t, j, k, m)
        c = 0.0 + 0.0im
        for i = 1:N
            c += ∂c_func(p, t, i, m) * aj_hi_aks[j, i, k]
        end
        adiff = a_diffs[j, k]
        return cis(t * adiff) * c
    end

    return ControllableSystem{
        Nothing,
        Nothing,
        typeof(new_c_func),
        typeof(new_∂c_func),
        Nothing,
    }(
        nothing,
        nothing,
        nothing,
        c_ls,
        c_ms,
        c_vs,
        new_c_func,
        new_∂c_func,
        nothing,
        nothing,
        false,
        d,
        true,
        sparse_tol,
    )
end

"""
    ControllableSystem(drift_op, basis_ops, RF_generator::Matrix, c_func, ∂c_func; <keyword arguments>)

Instantiate a ControllableSystem struct with a specified `RF_generator`. 

The `RF_generator` is a specification of the time-independent generator of the rotating
reference frame. i.e., if ``∂V(t)=A*V(t)`` where A is the generator of type ``Matrix``.

# Keyword Arguments
- `sparse_tol::Float`: A tolerance which defines a threshold to discard matrix elements
"""
function ControllableSystem(
    drift_op,
    basis_ops,
    RF_generator::Matrix,
    c_func,
    ∂c_func;
    sparse_tol = 1e-12,
)
    d = size(drift_op, 1)
    F = eigen(RF_generator)
    N = size(basis_ops, 1)
    c_ls = []
    c_ms = []
    c_vs = []
    as = F.values
    a_diffs = zeros(ComplexF64, d, d)
    aj_drift_aks = zeros(ComplexF64, d, d)
    aj_hi_aks = zeros(ComplexF64, d, N, d)
    for j = 1:d
        for k = 1:d
            aj = F.values[j]
            ak = F.values[k]
            a_diffs[j, k] = aj - ak
            aj_vec = @view F.vectors[:, j]
            ak_vec = @view F.vectors[:, k]
            aj_drift_ak = adjoint(aj_vec) * drift_op * ak_vec
            aj_drift_aks[j, k] = aj_drift_ak

            new_basis_op = sparse(aj_vec * adjoint(ak_vec))
            droptol!(new_basis_op, sparse_tol)
            for i = 1:N
                aj_hi_ak = adjoint(aj_vec) * basis_ops[i] * ak_vec
                aj_hi_aks[j, i, k] = aj_hi_ak
            end

            ls, ms, vs = findnz(new_basis_op)
            push!(c_ls, ls)
            push!(c_ms, ms)
            push!(c_vs, vs)
        end
    end

    function new_c_func(p, t, j, k)
        c = 0.0 + 0.0im
        diag_term = 0.0 + 0.0im
        if j == k
            diag_term = as[j]
        end
        for i = 1:N
            c += c_func(p, t, i) * aj_hi_aks[j, i, k]
        end
        adiff = a_diffs[j, k]
        aj_drift_ak = aj_drift_aks[j, k]
        return cis(t * adiff) * (aj_drift_ak + c + diag_term)
    end

    function new_∂c_func(p, t, j, k, m)
        c = 0.0 + 0.0im
        for i = 1:N
            c += ∂c_func(p, t, i, m) * aj_hi_aks[j, i, k]
        end
        adiff = a_diffs[j, k]
        return cis(t * adiff) * c
    end

    return ControllableSystem{
        Nothing,
        Nothing,
        typeof(new_c_func),
        typeof(new_∂c_func),
        Nothing,
    }(
        nothing,
        nothing,
        nothing,
        c_ls,
        c_ms,
        c_vs,
        new_c_func,
        new_∂c_func,
        nothing,
        nothing,
        false,
        d,
        true,
        sparse_tol,
    )
end

"""
    ControllableSystem(drift_op, basis_ops, RF_generator::Eigen, c_func, ∂c_func; <keyword arguments>)

Instantiate a ControllableSystem struct with a specified `RF_generator`. 

The `RF_generator` is a specification of the time-independent generator of the rotating
reference frame. i.e., if ``∂V(t)=A*V(t)`` where A is the generator of type ``LinearAlgebra.Diagonal``.

# Keyword Arguments
- `sparse_tol::Float`: A tolerance which defines a threshold to discard matrix elements
"""
function ControllableSystem(
    drift_op,
    basis_ops,
    RF_generator::Diagonal,
    c_func,
    ∂c_func,
)
    d_ls, d_ms, d_vs = findnz(drift_op)
    d = size(drift_op, 1)
    
    c_ls = [findnz(op)[1] for op in basis_ops]
    c_ms = [findnz(op)[2] for op in basis_ops]
    c_vs = [findnz(op)[3] for op in basis_ops]
    
    return ControllableSystem{
        typeof(d_ls),
        typeof(d_vs),
        typeof(c_func),
        typeof(∂c_func),
        typeof(RF_generator),
    }(
        d_ms, d_ls, d_vs,
        c_ls, c_ms, c_vs,
        c_func, ∂c_func,
        RF_generator,
        similar(RF_generator),
        true,
        d,
        false,
        0.0,
    )
end



"""
Instantiate a QOCProblem struct.
"""
struct QOCProblem
    Pc::Array{ComplexF64}
    Pc_dim::Int64
    Pa::Array{ComplexF64}
    Pa_dim::Int64
    target::Array{ComplexF64}
    control_time::Float64
end

"""
    QOCProblem(target, control_time, Pc, Pa)

Instantiate a QOCProblem struct with specified input operators and control time.

# Arguments
- `target::Array{ComplexF64}`: the target unitary operator in the computational subspace.
- `control_time::Float64`: The duration of the control.
- `Pc::Array{ComplexF64}`: A projector from the full unitary operator on the system to the computational subspace
- `Pa::Array{ComplexF64}`: A projector from the full unitary operator on the system to any ancillary subspace
"""
function QOCProblem(target, control_time, Pc, Pa)
    return QOCProblem(
        Pc,
        Int(round(tr(Pc))),
        Pa,
        Int(round(tr(Pa))),
        target,
        control_time
    )
end



"""
    g_sm(V, ∂U)

Computes infidelity via the squared modulus of the overlap between `V` and `U`.

This definition is perhaps the most standard used for quantum optimal control of unitaries.
"""
function g_sm(V, U; dim = 0, args...)
    if dim == 0
        dim = size(V, 1)
    end
    @assert(size(V) == size(U), "Unitaries are different sizes")
    return abs(1 - (1 / dim^2) * real(abs2(tr(adjoint(V) * U))))
end

"""
    ∂g_sm(V, ∂U)

Computes the partial derivative of `g_sm` w.r.t `U` evaluated at `∂U`.
"""
function ∂g_sm(V, U, ∂U; dim = 0, args...)
    if dim == 0
        dim = size(V, 1)
    end
    return (-2 / (dim^2)) * real(tr(adjoint(V) * ∂U) * tr(adjoint(U) * V))
end




"""
    GOAT_infidelity_reduce_map(sys, prob, goat_sol)

Maps the GOAT ODE solution to the objective function and gradient vector using an infidelity measure.

# Arguments
- `sys::ControllableSystem`: The controllable system.
- `prob::QOCProblem`: The QOCProblem
- `goat_sol::OrdinaryDiffEq.solution`: The solution to the GOAT equations of motion.
"""
function GOAT_infidelity_reduce_map(sys::ControllableSystem, prob::QOCProblem, goat_sol)
    d = sys.dim
    Pc = prob.Pc
    Pc_dim = prob.Pc_dim
    target = prob.target
    goatU = goat_sol.u[end]
    n_params = size(goatU, 1) ÷ d
    Ut = goatU[1:d, :]
    g = g_sm(Pc * target * Pc, Ut; dim = Pc_dim)

    ∂g_vec = Float64[]
    for i = 1:n_params-1
        ∂Ut = goatU[i*d+1:(i+1)*d, :]
        ∂g = ∂g_sm(Pc * target * Pc, Ut, ∂Ut; dim = Pc_dim)
        push!(∂g_vec, ∂g)
    end
    return g, ∂g_vec
end

"""
    SE_infidelity_reduce_map(sys, prob, SE_sol)

Maps Schrodinger ODE solution to the objective function using an infidelity measure.

# Arguments
- `sys::ControllableSystem`: The controllable system.
- `prob::QOCProblem`: The QOCProblem
- `SE_sol::OrdinaryDiffEq.solution`: The solution to the Schrodinger equation.
"""
function SE_infidelity_reduce_map(sys::ControllableSystem, prob::QOCProblem, SE_sol)
    d = sys.dim
    Pc = prob.Pc
    Pc_dim = prob.Pc_dim
    target = prob.target
    Ut = SE_sol.u[end]
    g = g_sm(Pc * target * Pc, Ut; dim = Pc_dim)
    return g
end


"""
Instantiate a QOCParameters struct.
"""
struct QOCParameters{A,B,C,D,E}
    ODE_options::NamedTuple
    SE_reduce_map::A
    GOAT_reduce_map::B
    optim_alg::C
    optim_options::D
    derivs_per_core::E
end

"""
    QOCParameters(ODE_options,SE_reduce_map, GOAT_reduce_map, optim_alg, optim_options; <keyword arguments> )

Instantiate a QOCParameters struct with specified parameters. 

# Arguments
- `ODE_options::NamedTuple`: The settings that will be input to the ODE solver.
- `SE_reduce_map::Function`: A function mapping the output from the Schrodinger equation to the objective value.
- `GOAT_reduce_map::Function`: A function mapping the output from the GOAT E.O.M.s to the objective value and gradient.
- `optim_alg::Optim Algorithm`: The specific optimization algorithm specified via Optim.jl
- `optim_options::Optim options`: The optimization algorithm options specified via Optim.jl
"""
function QOCParameters(
    ODE_options,
    SE_reduce_map,
    GOAT_reduce_map,
    optim_alg,
    optim_options;
    derivs_per_core = nothing,
)
    return QOCParameters{
        typeof(SE_reduce_map),
        typeof(GOAT_reduce_map),
        typeof(optim_alg),
        typeof(optim_options),
        typeof(derivs_per_core),
    }(
        ODE_options,
        SE_reduce_map,
        GOAT_reduce_map,
        optim_alg,
        optim_options,
        derivs_per_core,
    )
end



"""
    solve_GOAT_eoms(sys, opt_param_inds, Tmax, p; <keyword arguments>)

Integrate the Schrodinger equation for a specified time and control parameter set. 

# Arguments
- `sys::ControllableSystem`: The controllable system.
- `opt_param_inds::Vector{Int64}`: The vector of parameter indices specifying which gradients will be propogated.
- `Tmax::Float64`: The total contorl time.
- `p::Vector{Float64}`: The parameters which define the controlled evolution.
- `ODE_options`: The specification of the integrator settings from OrdinaryDiffEq.jl
"""
function solve_GOAT_eoms(
    sys::ControllableSystem,
    opt_param_inds::Vector{Int},
    Tmax::Float64,
    p::Vector{Float64};
    t0 = 0.0,
    args...,
)
    tspan = (t0, Tmax)
    g = make_GOAT_update_function(sys, opt_param_inds)
    u0 = make_GOAT_initial_state(sys.dim, opt_param_inds)
    prob = ODEProblem(g, u0, tspan, p)
    sol = solve(prob; args...)
    return sol
end




"""
    solve_GOAT_eoms_reduce(p, sys, prob, opt_param_inds, params::QOCParameters) 

Solves the GOAT eoms and outputs a objective function and gradient vector.

# Arguments
- `p`: The control parameter vector at which the objective and gradient is being calculated. 
- `sys::ControllableSystem`: The controllable system.
- `opt_param_inds`: The vector of parameter indices which determines which gradients are calculated.
- `param::QOCParameters`: The QOCParameters which provides the ODE_options.
"""
function solve_GOAT_eoms_reduce(
    p,
    sys::ControllableSystem,
    prob::QOCProblem,
    opt_param_inds,
    params::QOCParameters,
)
    T = prob.control_time
    goat_sol = solve_GOAT_eoms(sys, opt_param_inds, T, p; params.ODE_options...)
    out = params.GOAT_reduce_map(sys, prob, goat_sol)
    g = first(out)
    ∂gs = last(out)
    return g, ∂gs
end




"""
    parallel_GOAT_fg!(F, G, p, p_storage, opt_param_inds, sys, prob, params)

Parallelized computation of the objective and gradient for QOC with GOAT.

# Arguments
- `F`: The objective value
- `G`: The vector of gradients w.r.t. the control parameters `p`. 
- `p`: The control parameter vector.
- `p_storage`: A pre-allocated storage vector for current `p` values. 
- `opt_param_inds`: The vector of parameter indices which determines which gradients are calculated.
- `sys::ControllableSystem`: The controllable system.
- `prob::QOCProblem`: The QOCProblem
- `param::QOCParameters`: The QOCParameters.
"""
function parallel_GOAT_fg!(
    F,
    G,
    p,
    p_storage,
    opt_param_inds,
    sys::ControllableSystem,
    prob::QOCProblem,
    params::QOCParameters,
)
    T = prob.control_time
    p_storage[opt_param_inds] .= p # Update the storage vector with new parameters from optimization
    if G !== nothing
        num_params = size(p, 1)
        if params.derivs_per_core === nothing
            derivs_per_core = num_params
        else
            derivs_per_core = params.derivs_per_core
        end
        goat_param_indices =
            collect.(collect(Iterators.partition(opt_param_inds, derivs_per_core)))
        f = y -> solve_GOAT_eoms_reduce(p_storage, sys, prob, y, params)
        out = pmap(f, goat_param_indices)
        gs = first.(out)
        # @assert gs[1] ≈ gs[end] # Trivial sanity check
        for i = 1:size(goat_param_indices, 1)
            start = derivs_per_core * (i - 1) + 1
            stop = derivs_per_core * i
            ∂gs = last(out[i])
            G[start:stop] .= ∂gs
        end
        g = gs[1]
    else
        sol = solve_SE(sys, T, p_storage; params.ODE_options...)
        g = params.SE_reduce_map(sys, prob, sol)
    end

    if F !== nothing
        return g
    end

end



"""
    find_optimal_controls(p0, opt_param_inds, sys, prob, params)

Run the GOAT algorithm and find optimal controls.

# Arguments:
- `p0`: The initial guess of the optimal control parameters.
- `opt_param_inds`: Indices of p0 which specify which parameters to hold constant and which to optimize.
- `sys::ControllableSystem`: The controllable system.
- `prob::QOCProblem`: The quantum optimal control problem.
- `param::QOCParameters`: The quantum optimal control parameters.
"""
function find_optimal_controls(
    p0,
    opt_param_inds,
    sys::ControllableSystem,
    prob::QOCProblem,
    params::QOCParameters,
)
    p_storage = deepcopy(p0)
    fg!(F, G, x) = parallel_GOAT_fg!(F, G, x, p_storage, opt_param_inds, sys, prob, params)
    x0 = p0[opt_param_inds]
    res = Optim.optimize(Optim.only_fg!(fg!), x0, params.optim_alg, params.optim_options)
    return res
end







"""
    find_optimal_controls(p0, sys, prob, params)

Run the GOAT algorithm and find optimal controls.

# Arguments:
- `p0`: The initial guess of the optimal control parameters.
- `sys::ControllableSystem`: The controllable system.
- `prob::QOCProblem`: The quantum optimal control problem.
- `param::QOCParameters`: The quantum optimal control parameters.
"""
function find_optimal_controls(
    p0,
    sys::ControllableSystem,
    prob::QOCProblem,
    params::QOCParameters,
)
    fg!(F, G, x) = parallel_GOAT_fg!(F, G, x, sys, prob, params)
    res = Optim.optimize(Optim.only_fg!(fg!), p0, params.optim_alg, params.optim_options)
    return res
end


