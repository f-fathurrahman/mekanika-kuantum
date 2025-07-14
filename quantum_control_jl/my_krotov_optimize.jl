
function my_optimize_krotov(problem)

    @info "ENTER optimize_krotov"

    callback = get(problem.kwargs, :callback, (args...) -> nothing)
    if haskey(problem.kwargs, :update_hook) || haskey(problem.kwargs, :info_hook)
        msg = "The `update_hook` and `info_hook` arguments have been superseded by the `callback` argument"
        throw(ArgumentError(msg))
    end
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    # note: the default `check_convergence!` is a no-op. We still always check
    # for "Reached maximum number of iterations" in `update_result!`
    verbose = get(problem.kwargs, :verbose, false)
    skip_initial_forward_propagation =
        get(problem.kwargs, :skip_initial_forward_propagation, false)

    wrk = KrotovWrk(problem; verbose)

    ϵ⁽ⁱ⁾ = wrk.pulses0
    ϵ⁽ⁱ⁺¹⁾ = wrk.pulses1

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        println("Initial propagation")
        for (k, traj) in collect(enumerate(wrk.trajectories))
            println("Trajectory k = ", k)
            my_krotov_initial_fw_prop!(ϵ⁽ⁱ⁾, traj.initial_state, k, wrk)
        end
    end

    # TODO: if sigma, fw_storage0 = fw_storage
    update_result!(wrk, 0)
    info_tuple = callback(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    println("info_tuple = ", info_tuple)
    if !(isnothing(info_tuple) || isempty(info_tuple))
        push!(wrk.result.records, info_tuple)
    end

    i = wrk.result.iter  # = 0, unless continuing from previous optimization
    println("i = ", i)
    atexit_filename = get(problem.kwargs, :atexit_filename, nothing)
    println("atexit_filename = ", atexit_filename)
    # atexit_filename is undocumented on purpose: this is considered a feature
    # of @optimize_or_load
    if !isnothing(atexit_filename)
        set_atexit_save_optimization(atexit_filename, wrk.result)
        if !isinteractive()
            @info "Set callback to store result in $(relpath(atexit_filename)) on unexpected exit."
            # In interactive mode, `atexit` is very unlikely, and
            # `InterruptException` is handles via try/catch instead.
        end
    end
    # try-catch statement is removed
    while !wrk.result.converged
        i = i + 1
        println("Begin iteration = ", i)
        krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)
        update_result!(wrk, i)
        info_tuple = callback(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
        if !(isnothing(info_tuple) || isempty(info_tuple))
            push!(wrk.result.records, info_tuple)
        end
        check_convergence!(wrk.result)
        ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
    end


    finalize_result!(ϵ⁽ⁱ⁾, wrk)
    if !isnothing(atexit_filename)
        popfirst!(Base.atexit_hooks)
    end

    @info "EXIT my_optimize_krotov"

    return wrk.result

end


using MyQuantumPropagators.Storage: write_to_storage!

function my_krotov_initial_fw_prop!(ϵ⁽⁰⁾, ϕₖ, k, wrk)

    @info "ENTER my_krotov_initial_fw_prop"

    for propagator in wrk.fw_propagators
        propagator.parameters = IdDict(zip(wrk.controls, ϵ⁽⁰⁾))
    end
    #@infiltrate

    # set parameter?
    reinit_prop!(wrk.fw_propagators[k], ϕₖ; MyKrotov.transform_control_ranges)

    Φ₀ = wrk.fw_storage[k]
    (Φ₀ !== nothing) && write_to_storage!(Φ₀, 1, ϕₖ)
    # This will simply: Φ₀[1] <- ϕₖ


    N_T = length(wrk.result.tlist) - 1
    println("N_T = ", N_T)
    for n = 1:N_T
        Ψₖ = prop_step!(wrk.fw_propagators[k])
        # In case we have callback for propagation step
        if haskey(wrk.fw_prop_kwargs[k], :callback)
            local cb = wrk.fw_prop_kwargs[k][:callback]
            local observables = get(wrk.fw_prop_kwargs[k], :observables, _StoreState())
            cb(wrk.fw_propagators[k], observables)
        end
        (Φ₀ !== nothing) && write_to_storage!(Φ₀, n + 1, Ψₖ)
    end
end