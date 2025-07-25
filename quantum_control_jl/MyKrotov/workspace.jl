import MyQuantumControl
using MyQuantumControl.Controls: get_control_derivs
using MyQuantumControl.MyQuantumPropagators.Controls: get_controls, discretize_on_midpoints
using MyQuantumControl: init_prop_trajectory
using MyQuantumControl.MyQuantumPropagators.Storage: init_storage
using ConcreteStructs: @concrete

"""Krotov workspace.

The workspace is for internal use. However, it is also accessible in a
`callback` function. The callback may use or modify some of the following
attributes:

* `trajectories`: a copy of the trajectories defining the control problem
* `adjoint_trajectories`: The `trajectories` with the adjoint generator
* `kwargs`: The keyword arguments from the [`ControlProblem`](@ref) or the
  call to [`optimize`](@ref).
* `controls`: A tuple of the original controls (probably functions)
* `ga_a_int`: The current value of ``∫gₐ(t)dt`` for each control
* `update_shapes`: The update shapes ``S(t)`` for each pulse, discretized on
  the intervals of the time grid.
* `lambda_vals`: The current value of λₐ for each control
* `result`: The current result object
* `fw_storage`: The storage of states for the forward propagation
* `fw_propagators`: The propagators used for the forward propagation
* `bw_propagators`: The propagators used for the backward propagation
* `use_threads`: Flag indicating whether the propagations are performed in
  parallel.
"""
mutable struct KrotovWrk

    trajectories
    adjoint_trajectories
    kwargs
    controls
    # storage for controls discretized on intervals of tlist
    pulses0::Vector{Vector{Float64}}
    # second pulse storage: pulses0 and pulses1 alternate in storing the guess
    # pulses and optimized pulses in each iteration
    pulses1::Vector{Vector{Float64}}
    g_a_int::Vector{Float64}
    update_shapes::Vector{Vector{Float64}}
    lambda_vals::Vector{Float64}
    J_T_takes_tau::Bool  # Does J_T have a tau keyword arg?
    chi_takes_tau::Bool # Does chi have a tau keyword arg?
    # map of controls to options
    result

    #################################
    # Per trajectory:

    control_derivs
    fw_prop_kwargs::Vector{Dict{Symbol,Any}}
    bw_prop_kwargs::Vector{Dict{Symbol,Any}}
    fw_storage # forward storage array (per trajectory)
    fw_storage2 # forward storage array (per trajectory)
    bw_storage # backward storage array (per trajectory)
    fw_propagators
    bw_propagators
    use_threads::Bool

end


function KrotovWrk(problem::MyQuantumControl.ControlProblem; verbose=false)
    use_threads = get(problem.kwargs, :use_threads, false)
    trajectories = [traj for traj in problem.trajectories]
    N = length(trajectories)
    adjoint_trajectories = [adjoint(traj) for traj in problem.trajectories]
    controls = get_controls(trajectories)
    if length(controls) == 0
        error("no controls in trajectories: cannot optimize")
    end
    control_derivs = [get_control_derivs(traj.generator, controls) for traj in trajectories]
    tlist = problem.tlist
    kwargs = Dict(problem.kwargs)  # creates a shallow copy; ok to modify
    default_update_shape = get(problem.kwargs, :update_shape, t -> 1.0)
    default_lambda_a = convert(Float64, get(problem.kwargs, :lambda_a, 1.0))
    default_pulse_options = IdDict(
        c => Dict(:lambda_a => default_lambda_a, :update_shape => default_update_shape)
        for c ∈ controls
    )
    if :pulse_options in keys(kwargs)
        if :update_shape in keys(kwargs)
            @warn("`update_shape` is ignored due to given `pulse_options`")
        end
        if :lambda_a in keys(kwargs)
            @warn("`lambda_a=$(kwargs[:lambda_a])` is ignored due to given `pulse_options`")
        end
    else
        if (:update_shape ∉ keys(kwargs)) && (:lambda_a ∉ keys(kwargs))
            @warn "Using default pulse_options: (:lambda_a => 1.0, :update_shape => (t -> 1.0))"
        end
    end
    pulse_options = get(kwargs, :pulse_options, default_pulse_options)
    for c ∈ controls
        if c ∉ keys(pulse_options)
            error("pulse_options must be defined for all controls")
        end
    end
    update_shapes = [
        discretize_on_midpoints(pulse_options[control][:update_shape], tlist) for
        control in controls
    ]
    lambda_vals =
        [convert(Float64, pulse_options[control][:lambda_a]) for control in controls]
    if haskey(kwargs, :continue_from)
        @info "Continuing previous optimization"
        result = kwargs[:continue_from]
        if !(result isa KrotovResult)
            # account for continuing from a different optimization method
            result = convert(KrotovResult, result)
        end
        result.iter_stop = get(kwargs, :iter_stop, 5000)
        result.converged = false
        result.start_local_time = now()
        result.message = "in progress"
        pulses0 = [
            discretize_on_midpoints(control, tlist) for control in result.optimized_controls
        ]
    else
        result = KrotovResult(problem)
        pulses0 = [discretize_on_midpoints(control, tlist) for control in controls]
    end
    pulses1 = [copy(pulse) for pulse in pulses0]
    g_a_int = zeros(length(pulses0))
    # TODO: forward_storage only if g_b != 0
    fw_storage = [init_storage(traj.initial_state, tlist) for traj in trajectories]
    # TODO: second forward storage only if second order
    fw_storage2 = [init_storage(traj.initial_state, tlist) for traj in trajectories]
    bw_storage = [init_storage(traj.initial_state, tlist) for traj in trajectories]
    kwargs[:piecewise] = true  # only accept piecewise propagators
    _prefixes = ["prop_", "fw_prop_"]
    fw_prop_kwargs = [Dict{Symbol,Any}() for _ = 1:N]
    bw_prop_kwargs = [Dict{Symbol,Any}() for _ = 1:N]
    fw_propagators = [
        init_prop_trajectory(
            traj,
            tlist;
            verbose,
            _msg="Initializing fw-prop of trajectory $k",
            _prefixes,
            _filter_kwargs=true,
            _kwargs_dict=fw_prop_kwargs[k],
            kwargs...
        ) for (k, traj) in enumerate(trajectories)
    ]
    _prefixes = ["prop_", "bw_prop_"]
    bw_propagators = [
        init_prop_trajectory(
            traj,
            tlist;
            verbose,
            _msg="Initializing bw-prop of trajectory $k",
            _prefixes,
            _filter_kwargs=true,
            bw_prop_backward=true,  # will filter to `backward=true`
            _kwargs_dict=bw_prop_kwargs[k],
            kwargs...
        ) for (k, traj) in enumerate(adjoint_trajectories)
    ]
    J_T_takes_tau = false
    if haskey(kwargs, :J_T)
        J_T = kwargs[:J_T]
    else
        msg = "`optimize` for `method=Krotov` must be passed the functional `J_T`."
        throw(ArgumentError(msg))
    end
    J_T_takes_tau =
        hasmethod(J_T, Tuple{typeof(result.states),typeof(trajectories)}, (:tau,))
    if !haskey(kwargs, :chi)
        kwargs[:chi] = make_chi(J_T, trajectories)
    end
    chi = kwargs[:chi]
    chi_takes_tau =
        hasmethod(chi, Tuple{typeof(result.states),typeof(trajectories)}, (:tau,))
    return KrotovWrk(
        trajectories,
        adjoint_trajectories,
        kwargs,
        controls,
        pulses0,
        pulses1,
        g_a_int,
        update_shapes,
        lambda_vals,
        J_T_takes_tau,
        chi_takes_tau,
        result,
        control_derivs,
        fw_prop_kwargs,
        bw_prop_kwargs,
        fw_storage,
        fw_storage2,
        bw_storage,
        fw_propagators,
        bw_propagators,
        use_threads
    )
end
