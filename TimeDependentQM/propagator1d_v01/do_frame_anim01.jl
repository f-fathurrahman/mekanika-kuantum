using LinearAlgebra
using Printf

using Plots
using PlotThemes
theme(:dark)

include("../build_D2_matrix_3pt.jl")
include("../build_D2_matrix_11pt.jl")

# A free Hamiltonian (no potential)
function build_Ham(N::Int64, Δx::Float64)
    D2 = build_D2_matrix_11pt(N, Δx)
    return -0.5*D2
end

# Using Ha unit: ħ = 1, m = 1

# Analytical propagation of Gaussian wavepacket
function gaussian_wavepacket(x, t, σ, k0; x0=0.0)
    Ω = σ^2
    v = k0
    oneo = (1 + im*Ω*t)
    pw_f = exp( im*k0*(x - v*t - x0) ) # shift by x0 ?
    gauss_f = exp( -σ^2/2 * (x - v*t - x0)^2 / oneo )
    pre_f = sqrt(σ/( sqrt(π)*oneo) )
    f = pre_f * gauss_f * pw_f
    return f
end

function main()

    a = -25.0
    b =  25.0

    σ = sqrt(0.1) # higher value will give sharper peak
    k0 = 5.0  # higher value will move the wavepacket faster
    x0 = -10.0

    k_grid_min = 2*k0 * 10
    Δx = 2π/k_grid_min
    println("Δx min = ", Δx)
    Nx = floor(Int64, (b - a)/Δx) + 2
    println("Calculated Nx = ", Nx)

    if Nx < 100
        Nx = 200
        println("Nx is set to ", Nx)
    end

    # NOTE: This is different from the one used in Varga
    x = zeros(Float64, Nx)
    Δx = (b - a)/(Nx - 1)
    for i in 1:Nx
        x[i] = a + (i-1)*Δx # calculate x grid point here
    end

    t_last = 4.5
    Nt = 200
    Δt = t_last/(Nt - 1)
    
    println("Δt = ", Δt)
    println("k0*t_last = ", k0*t_last)
    println("2π/k0 = ", 2π/k0)
    println("Δx    = ", Δx)
    println()
    println("k0     = ", k0)
    println("k grid = ", 2π/Δx)
    println("CFL = ", k0*Δt/Δx)

    # Calculate the Hamiltonian
    H = build_Ham(Nx, Δx)
    # The propagator, we are using inverse directly here
    U = inv(I + im*Δt/2 * H)*(I - im*Δt/2 * H)

    # initial wave function
    t = 0.0
    ψ = zeros(ComplexF64, Nx)
    for i in 1:Nx
        ψ[i] = gaussian_wavepacket(x[i], t, σ, k0, x0=x0)
    end

    # Check norm
    @printf("Initial norm of the wavepacket: %18.10f\n", dot(ψ, ψ)*Δx)

    plot(x, abs.(ψ), label="Abs", linewidth=2, dpi=150)
    title!("Initial state")
    ylims!(0.0, 0.8)
    xlims!(a, b)
    savefig("IMG_00000.png")


    plot(x, real.(ψ), label="Abs", linewidth=2, dpi=150)
    title!("Initial state")
    ylims!(0.0, 0.8)
    savefig("IMG_00000_real.png")

    plot(x, imag.(ψ), label="Abs", linewidth=2, dpi=150)
    title!("Initial state")
    ylims!(0.0, 0.8)
    savefig("IMG_00000_imag.png")


    t = 0.0
    ψ_ana = zeros(ComplexF64, Nx)
    ψ_new = zeros(ComplexF64, Nx)
    for itime in 1:Nt
        t += Δt
        @views ψ_new[:] = U * ψ[:] # numerical
        # analytic
        for i in 1:Nx
            ψ_ana[i] = gaussian_wavepacket(x[i], t, σ, k0, x0=x0)
        end
        ψ[:] = ψ_new[:]
        #
        plot(x, abs.(ψ_new), label="Numerical", linewidth=2, dpi=150)
        plot!(x, abs.(ψ_ana), label="Analytical", linewidth=2, dpi=150)
        title!("Final state")
        ylims!(0.0, 0.8)
        xlims!(a, b)
        filename = @sprintf("IMG_%05d.png", itime)
        savefig(filename)
    end
    println("Final t = ", t)

    @printf("After numerical propagation: norm of the wavepacket: %18.10f\n", dot(ψ,ψ)*Δx)

end
main()
