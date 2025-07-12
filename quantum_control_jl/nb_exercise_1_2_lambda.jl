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
# # Population Transfer in a Three-Level-System with STIRAP

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
# This notebook introduces the so-called Lambda-system, a three-level system which serves as a slightly more advanced version of the two-level system discussed in [Exercise I.1](jl_exercise_1_1_TLS.ipynb). In a Lambda system the levels 1 and 2 and the levels 2 and 3 are coupled, but the levels 1 and 3 are not. Drawing the level schemes with the corresponding couplings then looks a lot like the letter $\Lambda$, which you will see in the course of this notebook! We will consider two pulses interacting with a lambda system, one for each of the two couplings. While a transition from levels 1 to 3 is not directly allowed, it is possible via the intermediate level 2. At the end of this notebook you will learn, that it is even possible to transfer population from level 1 to level 3 without ever populating the intermediate level 2!

# %%
Pkg.activate("quantumcontrol", shared=true)

# %% [markdown]
# ## Setup

# %% [markdown]
# We load some of the Julia packages we will use in the following

# %%
using LinearAlgebra
using QuantumPropagators: hamiltonian, propagate
using OrdinaryDiffEq

# %% [markdown]
# We'll also set some defaults for `Plots`, like increasing the default line width for better readability.

# %%
using Plots, PlotThemes

# %%
theme(:dark)

# %%
# Some utilities for showing hints and solutions
include(joinpath("tutorial_utils", "exercise_1_lambda.jl"));

# %% [markdown]
# ## The STIRAP Hamiltonian

# %% [markdown]
# We consider a three level system with levels labeled 1,2,3. The levels 1-2 and 2-3 are connected via the corresponding transition dipole moments $d_{1,2}$ and $d_{2,3}$. An electric dipole transition between levels 1 and 3 is forbidden, leading to $d_{1,3}=0$.
#
# The system is initialized in the state $\ket{1}$. The goal is to transfer the population to the state $\ket{3}$ by means of coherent control. Furthermore, we choose as an additional goal to avoid populating level 2 as much as possible. This could be physically motivated for example by the state $\ket{2}$ being subject to unwanted decay.
#
# We will use two light fields as controls, targeting the transition 1-2 and 2-3. These are usually called the *Pump* and *Stokes* pulse, respectively.
#
# The Hamiltonian in the lab frame takes the form
#
# \begin{align}\op{H}_{\text{lab}} &= \begin{pmatrix}\varepsilon_{1} & d^{*}_{1,2}E(t)\\ d_{1,2}E(t) & \varepsilon_{2} & d^{*}_{2,3}E(t)\\ & d_{2,3}E(t) & \varepsilon_{3} \end{pmatrix},\end{align}
#
# where $E(t)=E_{\text{P}}(t) + E_{\text{S}}(t)$ is the time-dependent electric control field which consists of both the pump and stokes fields $E_{\text{P/S}}(t)$. The fields have frequencies $\omega_{\text{P}}$ resp. $\omega_{\text{S}}$ and take the form
# \begin{align} E_{P/S}(t)&=S_{\text{P/S}}(t)\cos \left(\omega_{\text{P/S}}t\right), \end{align}
#
# with $S_{\text{P/S}}(t)$ the *shape functions*. We go to the rotating frame defined by the following unitary matrix
# \begin{align} \op{U} &= e^{i\epsilon_{1}t}\begin{pmatrix}1\\ & e^{i\omega_{P}t}\\ & & e^{i\left(\omega_{P}-\omega_{S}\right)t} \end{pmatrix}. \end{align}
#
# Under this transformation the Hamiltonian transforms as
# \begin{align} \op{H}_{\text{rot}}=\op{U}\op{H}_{\text{lab}}\op{U}^{\dagger}+i\dot{\op{U}}\op{U}^{\dagger}. \end{align}
#
# After performing a rotating wave approximation, in which fast oscillating terms like $e^{i(\omega_{\text{P}} \pm \omega_{\text{S}})t}$ are neglected, we obtain the Hamiltonian most commonly found in literature as the starting point of the STIRAP protocol,
#
# \begin{align} \op{H}_{\text{STIRAP}}=\begin{pmatrix}0 & \frac{1}{2}\Omega^{*}_{\text{P}}(t)\\
# \frac{1}{2}\Omega_{\text{P}}(t) & \Delta_\text{P} & \frac{1}{2}\Omega^{*}_{\text{S}}(t)\\
# & \frac{1}{2}\Omega_{\text{S}}(t) & \Delta_\text{P}-\Delta_\text{S} \end{pmatrix},
# \end{align}
#
# where $\Delta_{\text{P}} \equiv (\varepsilon_2 - \varepsilon_1) - \omega_\text{P}$ and $\Delta_S \equiv (\varepsilon_2 - \varepsilon_3) - \omega_\text{S}$ are called *one-photon detunings* and $\delta \equiv \Delta_\text{P} - \Delta_\text{S}$ is the *two-photon detuning*. We also defined $\Omega_{\text{P}}(t) = d_{1,2} S_\text{P}(t)$ resp. $\Omega_{\text{S}}(t) = d_{2,3} S_\text{S}(t)$, the time-dependent *Rabi frequencies*.
#
# The level scheme, pulses and detunings are illustrated in the following figure:
#
# <center><img src="../figures/lambda_system_levels.png" alt="Lambda system considered in this notebook" width="500"></center>
#
# As mentioned in the introduction, this configuration is called Lambda system because the level scheme resembles the letter $\Lambda$!

# %% [markdown]
# ## The dark state and STIRAP protocol
#
# In the following, we assume the controls to be two-photon resonant, i.e. $\delta = \Delta_\text{P} - \Delta_\text{S} = 0$. We write $\Delta \equiv \Delta_\text{P} = \Delta_\text{S}$.
#
# The three instantaneous eigenstates of $\hat{H}_\text{STIRAP}$ can then be expressed as
#
# \begin{align}
# \ket{\varphi_{+}} & =\sin\varphi\sin\theta\ket{1}+\cos\varphi\ket{2}+\sin\varphi\cos\theta\ket{3}\\
# \ket{\varphi_{-}} & =\cos\varphi\sin\theta\ket{1}-\sin\varphi\ket{2}+\cos\varphi\cos\theta\ket{3}\\
# \ket{\varphi_{0}} & =\cos\theta\ket{1}-\sin\theta\ket{3}
# \end{align}
#
# with corresponding eigenvalues $\varepsilon_\pm(t) = \frac{\Delta(t)}{2}\pm\frac{1}{2}\sqrt{\Delta^2(t) + \Omega_0^2(t)}$ and $\Omega_0(t)=\sqrt{\Omega_\text{P}^2(t) + \Omega_\text{S}^2(t)}$. The time-dependent angles are given via
#
# \begin{align}
# \tan \theta(t) &= \frac{\Omega_\text{P}(t)}{\Omega_\text{S}(t)} \\
# \tan 2\varphi(t) &= \frac{\Omega_0 (t)}{\Delta(t)}
# \end{align}
#
# Notice that the eigenstate $\ket{\varphi_0}$ does *not* contain any amplitude in state $\ket{2}$. It is therefore immune to any decay affecting $\ket{2}$. Since it is not directly involved in the process, it is called a *dark state*.
#
# Since the dark state is a superposition of $\ket{1}$ and $\ket{3}$ where the respective amplitudes depend on the control amplitudes, we can use the dark state to transfer population from 1 to 3 by having the system adiabatically follow the dark state.
#
# Initially (at $t=-\infty$), we want the system to be in state $\ket{1}$. For $\ket{\phi_0}$, this corresponds to $\theta = 0$ or $\Omega_\text{P} / \Omega_\text{S} \rightarrow 0$.
# Eventually (at $t=\infty$), we want the system to be in state $\ket{3}$, demanding $\theta = \pi/2$ or $\Omega_\text{S} / \Omega_\text{P} \rightarrow 0$.
#
# This indicates a counter-intuitive ordering of the two pulses: At first we need to apply the Stokes pulse targeting the 2-3 transition and thereafter we need to apply the Pump pulse on the 1-2 transition.
#
# On top of the correct order of the pulses, we need to make sure that the protocol is slow enough such the population always follows the instantaneous eigenstates. This is known as *adiabaticity*. One can show [1] that the system follows its instantaneous eigenstates adiabatically if
#
# \begin{align}
# \Omega_{0}(t)\gg & \left|\dot{\theta}(t)\right|\\
# \Leftrightarrow1\gg & \left|\frac{\Omega_{\text{S}}(t)\dot{\Omega}_{\text{P}}(t)-\Omega_{P}(t)\dot{\Omega}_{S}(t)}{\Omega_{0}^{3}(t)}\right|
# \end{align}
#
# The above condition can be interpreted as follows: The duration of the pulses can be estimated as $T \approx 1/ |\dot{\theta}|$. The adiabaticity condition can then be written as $\Omega_0 T \gg 1$. Note that $\Omega_0$ is essentially the overall strength of the two pulses. In short: The protocol follows the dark state for very long and very strong pulses.
#
# [1] http://arxiv.org/abs/1605.00224v2

# %% [markdown]
# ## Implementation using Gaussian pulses
#
# Now let's consider a numerical implementation of the previous theoretical findings where we model the two pulses with Gaussian envelopes. First, we define the pulse shape functions and another function returning the STIRAP Hamiltonian.
#
# We also define the time interval `t_list` for the calculation and its initial state `psi_0`.

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

* `tₚ`: center of the Gaussian used for the pump pulse
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
function H_stirap(Ωₚ, Ωₛ; Δ, d₁₂, d₂₃)

    H₀ = zeros(3,3)
    H₀[2,2]=Δ

    H₁₂ = zeros(3,3)
    H₁₂[1,2] = d₁₂ / 2
    H₁₂ += H₁₂'

    H₂₃ = zeros(3,3)
    H₂₃[2,3] = d₂₃ / 2
    H₂₃ += H₂₃'

    return hamiltonian(H₀, (H₁₂, Ωₚ), (H₂₃, Ωₛ))
end

tlist = collect(range(-250., 250.; length=401))

psi_0 = ComplexF64[1, 0, 0]

# %% [markdown]
# ## Problem 1: Choosing parameters for STIRAP

# %% [markdown]
# Now we can specify the parameters for the calculation. Below you will find a
# set of parameters that does not result in a perfect STIRAP protocol. Play
# around with the parameters and try to create the ideal dynamics described above. In fact, the initial parameters we provided for you are already almost correct. There is just one small issue - can you find it?

# %%
Δ = 1.    # Single-photon detuning
d₁₂ = 25. # Electric-dipole moment for 1-2 transition / Pump Strength
d₂₃ = 25. # Electric-dipole moment for 2-3 transition / Stokes Strength

tₚ = -50.  # Center of the Pump pulse
tₛ = 51. # Center of the Stokes pulse
σₚ = 50.  # Width of the Pump pulse
σₛ = 50.  # Width of the Stokes pulse

Ωₚ(t) = pump_shape(t; tₚ, σₚ)
Ωₛ(t) = stokes_shape(t; tₛ, σₛ)

H = H_stirap(Ωₚ, Ωₛ; Δ, d₁₂, d₂₃);

# %%
states = propagate(psi_0, H, tlist; method=OrdinaryDiffEq, storage=true)

# %%
using QuantumPropagators.Controls: discretize

function plot_population_and_pulses(tlist, states, ϵₚ, ϵₛ)
    Eₚ = discretize(ϵₚ, tlist)
    Eₛ = discretize(ϵₛ, tlist)
    max_E = maximum(vcat(Eₚ, Eₛ))
    pul_plt = plot(tlist, abs.([Eₚ Eₛ]) ./ max_E; label=["|ϵₚ|" "|ϵₛ|"], color=["#1f77b4" "#ff7f0e"], xlabel="time", ylabel="pulse")
    stt_plt = plot(tlist, abs2.(states)'; label=["|1⟩" "|2⟩" "|3⟩"], ls=[:solid :solid :solid], color=["#1f77b4" "#ff7f0e" "#2ca02c"], xlabel="time", ylabel="population")
    fig = plot(pul_plt, stt_plt, size=(800,400))
    infidelity = 1 - abs2.(states[3,end])
    return fig
end

plot_population_and_pulses(tlist, states, Ωₚ, Ωₛ)

# %%
problem_1.solution

# %% [markdown]
# ## Next Steps

# %% [markdown]
# After learning about the basics of the STIRAP protocol, you can now try out to find the optimal pulse parameters numerically with optimal control! [Example II.2](jl_exercise_2_2_lambda.ipynb) explains how to find the proper parameters to achieve the population inversion discussed above with a gradient-free optimization. [Example III.2](jl_exercise_3_2_lambda.ipynb) does the same with a gradient-based approach - Krotov's method.
