# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.11.6
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# # Parameter Optimization for STIRAP

# %% [markdown]
# $\newcommand{tr}[0]{\operatorname{tr}}
# \newcommand{diag}[0]{\operatorname{diag}}
# \newcommand{abs}[0]{\operatorname{abs}}
# \newcommand{pop}[0]{\operatorname{pop}}
# \newcommand{aux}[0]{\text{aux}}
# \newcommand{opt}[0]{\text{opt}}
# \newcommand{tgt}[0]{\text{tgt}}
# \newcommand{init}[0]{\text{init}}
# \newcommand{lab}[0]{\text{lab}}
# \newcommand{rwa}[0]{\text{rwa}}
# \newcommand{bra}[1]{\langle#1\vert}
# \newcommand{ket}[1]{\vert#1\rangle}
# \newcommand{Bra}[1]{\left\langle#1\right\vert}
# \newcommand{Ket}[1]{\left\vert#1\right\rangle}
# \newcommand{Braket}[2]{\left\langle #1\vphantom{#2} \mid
# #2\vphantom{#1}\right\rangle}
# \newcommand{op}[1]{\hat{#1}}
# \newcommand{Op}[1]{\hat{#1}}
# \newcommand{dd}[0]{\,\text{d}}
# \newcommand{Liouville}[0]{\mathcal{L}}
# \newcommand{DynMap}[0]{\mathcal{E}}
# \newcommand{identity}[0]{\mathbf{1}}
# \newcommand{Norm}[1]{\lVert#1\rVert}
# \newcommand{Abs}[1]{\left\vert#1\right\vert}
# \newcommand{avg}[1]{\langle#1\rangle}
# \newcommand{Avg}[1]{\left\langle#1\right\rangle}
# \newcommand{AbsSq}[1]{\left\vert#1\right\vert^2}
# \newcommand{Re}[0]{\operatorname{Re}}
# \newcommand{Im}[0]{\operatorname{Im}}$

# %% [markdown]
# In this notebook we return to the three-level-system in a Lambda configuration introduced in [Exercise I.2](jl_exercise_1_2_lambda.ipynb). Our goal remains to achieve population transform from level 1 to level 3 without populating the intermediate level 2. In this notebook you will learn how to use gradient-free parameter optimization for this purpose. A particularly important part of the optimization is the definition of an appropriate optimization functional which incorporates the goal to avoid populating level 2 as an additional condition.

# %% [markdown]
# ## Setup

# %%
Pkg.activate("quantumcontrol", shared=true)

# %%
using LinearAlgebra
using QuantumPropagators: hamiltonian, propagate, discretize
using OrdinaryDiffEq
using Optimization, OptimizationNLopt
using ComponentArrays: ComponentVector
using UnPack: @unpack

# %%
# Some utilities for showing hints and solutions
include(joinpath("tutorial_utils", "exercise_2_lambda.jl"));

# %% [markdown]
# We'll again set some defaults for `Plots`, like increasing the default line width for better readability.

# %%
using Plots, PlotThemes

# %%
theme(:dark)

# %% [markdown]
# ## Model

# %% [markdown]
# As a reminder, the STIRAP Hamiltonian is given by the following expression
#
# \begin{align} H_{\text{STIRAP}}=\begin{pmatrix}0 & \frac{1}{2}\Omega^{*}_{\text{P}}(t)\\
# \frac{1}{2}\Omega_{\text{P}}(t) & \Delta_\text{P} & \frac{1}{2}\Omega^{*}_{\text{S}}(t)\\
# & \frac{1}{2}\Omega_{\text{S}}(t) & \Delta_\text{P}-\Delta_\text{S} \end{pmatrix}
# \end{align}
#
# and gives rise to a level scheme in a Lambda configuration
#
# <center><img src="../figures/lambda_system_levels.png" alt="Lambda system considered in this notebook" width="500"></center>
#
# We now use parameter optimization to find the right pulse shapes, assuming once again a Gaussian pulse shape for simplicity. Furthermore we assume to be on two-photon resonance and fix the single-photon detuning at $\Delta \equiv 
# \Delta_\text{P} = 1$ which defines a reference energy. Therefore there are 3 parameters for each pulse which can be individually tuned: Its temporal position, its temporal width and its strength resulting. All in all, we have a total of 6 parameters.
#
# Our functional will thus take a list of these 6 parameters as an input to calculate how close we are to the following two goals:
#
#  1. At final time, all population should be in state 3.
#  2. Throughout the evolution, population in 2 should be kept as close to zero as possible.
#
# A straightforward approach to construct functionals which take multiple physical goals into account is to simply sum up functionals for each of the individual goals. In our case this is achieved with the following definition,
#
# \begin{align}
#     \mathcal{F} &= \Braket{3}{\psi (t_f)} - \frac{1}{T} \int_{t_i}^{t_f} \Braket{2}{\psi(t)} dt
# \end{align}
#
# where $t_i$, $t_f$ is the initial resp. final time and $T=t_f - t_i$ is the
# total duration of the protocol. This fidelity has a maximum value of $1$ corresponding to both goals being achieved perfectly, i.e., the final state of the evolution is $\ket{3}$ and there is no poulation in $\ket{2}$ during the
# entire protocol. Since we frame our optimizations as minimizations, we thus attempt to minimize the functional $1 - \mathcal{F}$ in the following.

# %% [markdown]
# ## Shape functions

# %% [markdown]
# First we define the pulse shape functions and the Hamiltonian.

# %%
"""
Gaussian shape function centered around t₀ with width σ

# Keyword Arguments

* `t₀`: center of the Gaussian
* `σ`: width of the Gaussian
"""
function gaussian_shape(t; t₀, σ)
    return @. exp(-0.5*(t - t₀)^2/σ^2) / √(2π*σ)
end

"""
Shape function for the pump pulse

# Keyword Arguments

* `tₚ`: center of the Gaussian used fro the pump pulse
* `σₚ`: width of the Gaussian
"""
function pump_shape(t; tₚ, σₚ)
    return gaussian_shape(t; t₀=tₚ, σ=σₚ)
end

"""
Shape function for the Stokes pulse

# Keyword Arguments

* `tₛ`: center of the Gaussian used for the Stokes pulse
* `σₛ`: width of the Gaussian
"""
function stokes_shape(t; tₛ, σₛ)
    return gaussian_shape(t; t₀=tₛ, σ=σₛ)
end

"""
Function returning the STIRAP Hamiltonian
"""
function H_stirap(Ωₚ, Ωₛ; Δ)

    H₀ = zeros(3,3)
    H₀[2,2]=Δ

    H₁₂ = zeros(3,3)
    H₁₂[1,2] = 1 / 2
    H₁₂ += H₁₂'

    H₂₃ = zeros(3,3)
    H₂₃[2,3] = 1 / 2
    H₂₃ += H₂₃'

    return hamiltonian(H₀, (H₁₂, Ωₚ), (H₂₃, Ωₛ))
end

tlist = collect(range(-250., 250.; length=101))

psi_0 = ComplexF64[1, 0, 0];

# %% [markdown]
# We also define a function which simulates the dynamics under a given pulse (e.g. the guess or optimized field) to visualize the pulses and the dynamics:

# %%
function evolve_and_plot_parameterized_pulse(x; Δ=1.0, T=250, nt=101, psi_0=psi_0)
    if isempty(x)
        error("It seems you are still using the (empty) guess. Make sure to fill it with the guess for the paramters `tₚ`, `tₛ`, `σₚ`, `σₛ`, `Ωₚ⁰`, `Ωₛ⁰`.!")
    end

    tₚ, tₛ, σₚ, σₛ, Ωₚ⁰, Ωₛ⁰ = x
    Ωₚ(t) = Ωₚ⁰ * pump_shape(t; tₚ=tₚ, σₚ=σₚ)
    Ωₛ(t) = Ωₛ⁰ * stokes_shape(t; tₛ=tₛ, σₛ=σₛ)
    tlist = collect(range(-T, T, nt))
    H = H_stirap(Ωₚ, Ωₛ; Δ=Δ);
    Ψs = propagate(psi_0, H, tlist; method=OrdinaryDiffEq, storage=true)

    Eₚ = discretize(Ωₚ, tlist)
    Eₛ = discretize(Ωₛ, tlist)
    max_E = maximum(vcat(Eₚ, Eₛ))
    pul_plt = plot(tlist, abs.([Eₚ Eₛ]) ./ max_E; label=["|Ωₚ|" "|Ωₛ|"], color=["#1f77b4" "#ff7f0e"], xlabel="time", ylabel="pulse")
    stt_plt = plot(tlist, abs2.(Ψs)'; label=["|0⟩" "|1⟩" "|2⟩"], ls=[:solid :solid :solid], color=["#1f77b4" "#ff7f0e" "#2ca02c"], xlabel="time", ylabel="population")
    plot(pul_plt, stt_plt, size=(800,400))
end

# %% [markdown]
# ## Parameter optimization

# %% [markdown]
# In the following, we implement the functional introduced above and optimize it via the `NLopt` package.
#
# The input for our functional are the two blackman pulses, which are parametrized by three values each: position of the maximum $t_{P/S}$, amplitude $\Omega^{(0)}_{P/S}$, and width $\sigma_{P/S}$.

# %%
function J_loss(x, constants)
    tₚ, tₛ, σₚ, σₛ, Ωₚ⁰, Ωₛ⁰ = x
    @unpack Δ, T, nt, psi_0 = constants
    Ωₚ(t) = Ωₚ⁰ * pump_shape(t; tₚ=tₚ, σₚ=σₚ)
    Ωₛ(t) = Ωₛ⁰ * stokes_shape(t; tₛ=tₛ, σₛ=σₛ)

    tlist = collect(range(-T, T; length=nt))
    H = H_stirap(Ωₚ, Ωₛ; Δ=Δ);
    Ψs = propagate(psi_0, H, tlist; method=OrdinaryDiffEq, storage=true)

    ovlp3 = abs(Ψs[3,end]) # term ⟨3|ψ(T)⟩
    int_ket2 = sum(abs.(Ψs[2,:]))/size(Ψs,2) # term ∝ ∫ ⟨2|ψ(T)⟩ dt

    return 1 - ovlp3 + int_ket2
end

# %% [markdown]
# ### Problem 1: Choosing a guess pulse

# %% [markdown]
# In order for the optimization to be succesful, we need to specify bounds for
# the optimization parameters. Moreover, the result of the optimization heavily
# depends on the initial guess parameters. Try finding a set of initial
# parameters that will converge to 99% fidelity!
#
# <!-- A good guess is [10,-50, 50,50, 25,25] -->

# %% [markdown]
# We begin by defining a set of guess parameters and lower/upper bounds for the optimization:

# %%
problem_1.solution

# %%
bounds_lower = [-100, -100, 10, 10, 0, 0]
bounds_upper = [100, 100, 80, 80, 60, 60]

guess = [30, -30, 50, 50, 25, 25];

# %%
prob = OptimizationProblem(
    J_loss,
    guess,
    (; Δ=1.0, T=250, nt=101, psi_0=psi_0);  # this is a NamedTuple, which will be accessible via prob.p
    lb=bounds_lower,
    ub=bounds_upper,
    stopval=(1-0.995), # below which error to stop the optimization
);

# %% [markdown]
# We can check the quality of the guess pulse:

# %%
evolve_and_plot_parameterized_pulse(guess; prob.p...)

# %% [markdown]
# Finally, we can check the value of the functional:

# %%
J_loss(guess, prob.p)

# %% [markdown]
# ### Optimization

# %% [markdown]
# In the optimization, we want to keep track of the fidelity after each iteration. To this end, we define a ["callback" function](https://docs.sciml.ai/Optimization/stable/API/solve/#CommonSolve.solve-Tuple%7BOptimizationProblem,%20Any%7D) that the optimizer will run after each step.

# %%
n = 0
function callback(state, loss_val)
    global n
    fid = 1 - loss_val
    n += 1
    print("Iteration: $(n), current fidelity $(round(fid; digits=6))\r")
    return false
end

# %% [markdown]
# Lastly, we call `Optimization.solve` to run the optimization:

# %%
obtained_fidelities = Float64[];
res = Optimization.solve(prob, NLopt.LN_NELDERMEAD(); maxiters=1000, callback)
println("\n\nHighest fidelity reached: $(round((1 - res.objective) * 100; digits=1))%")
if res.objective < 5.0e-3
    println("\tcongratulations, you have obtained population inversion!")
else
    println("\tbad guess, please try again!")
end

# %%
# problem_1.hint

# %%
# problem_1.solution

# %% [markdown]
# ## Analyze optimization results

# %% [markdown]
# Finally, let's have a look at the optimized results:

# %%
evolve_and_plot_parameterized_pulse(res.u; Δ=1.0, T=250, nt=1001, psi_0=psi_0)

# %% [markdown]
# ## Next steps

# %% [markdown]
# To go another step up in system complexity with parameter optimization we recommend [Exercise II.3](py_exercise_2_3_chiral.ipynb) which discusses the very interesting case of three-wave-mixing in a chiral molecule modelled by a three-level system. Alternatively, if you are interested in optimization with a gradient-based approach, we recommend to have a look at [Exercise III.2](py_exercise_3_2_lambda.ipynb)  in which Krotov's method is used for the opimization you studied in this notebook.
