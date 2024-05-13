using LinearAlgebra

using Plots, PlotThemes
theme(:dark)

include("../D2_matrices.jl")

function gaussian_wavepacket(x; x0=0.0, k=5.0, σ=0.1)
    g = sqrt(1 / sqrt(π) / σ) * exp(-(x - x0)^2 / 2 / σ^2)
    return exp(im * k * (x - x0) ) * g
end


function propagate_CrankNicolson(
    psi0, Vx, xgrid, Δt; Nt=100, print_norm=false
)
    # Crank-Nicolson method for the 1D Schrodinger equation.
    # No. of spatial grid points
    Nx  = size(xgrid, 1)
    Δx = xgrid[2] - xgrid[1]
    # for periodic grid the grid is slightly different

    # The Hamiltonian matrix
    H = -0.5*build_D2_matrix_11pt(Nx, Δx) + diagm(Vx)
    
    # the two unitary matrices
    # The propagator, we are using inverse directly here
    U = inv(I + im*Δt/2 * H)*(I - im*Δt/2 * H)

    psi_t = zeros(ComplexF64, Nx, Nt)
    psi_t[:,1] = psi0[:]

    for itime in 1:Nt-1
        psi_t[:,itime+1] = U*psi_t[:,itime]
    end
    return psi_t
end


#function main()
    xmin = -10.0
    xmax =  10.0
    Nx = 999
    xgrid = range(xmin, stop=xmax, length=Nx)
    dx = xgrid[2] - xgrid[1]

    # the gaussian wavepacket as initial wavefunction
    x0 = -5.0
    k0 = 5.0
    σ_x = 1.0  # larger value give wider wave
    psi0 = zeros(ComplexF64, Nx)
    for i in 1:Nx
        psi0[i] = gaussian_wavepacket(xgrid[i], x0=x0, k=k0, σ=σ_x)
    end

    # Potential
    Vx = zeros(Float64, Nx)

    total_time = 2.0
    # time step in atomic units, 1 a.u. = 24.188 as
    dt = 2 / (2 * π / σ_x + k0)^2
    dt0 = 0.001
    println("Time step: $(dt)")
    println("Actual time step: $(dt0)")
    # No. of temporal grid points
    Nt0 = round(Int64, total_time / dt0) + 1
    println("Nt0 = $(Nt0)")
    psi_t = propagate_CrankNicolson(psi0, Vx, xgrid, dt0, Nt=Nt0)

#end
#main()
