using Plots
using LinearAlgebra

function plot_bloch_state(θ::Real, ϕ::Real; title="")
    # Convert spherical coordinates to Cartesian
    x = sin(θ) * cos(ϕ)
    y = sin(θ) * sin(ϕ)
    z = cos(θ)
    
    # Create sphere
    u = range(0, 2π, length=50)
    v = range(0, π, length=50)
    xs = [cos(ϕ)*sin(θ) for θ in v, ϕ in u]
    ys = [sin(ϕ)*sin(θ) for θ in v, ϕ in u]
    zs = [cos(θ) for θ in v, ϕ in u]
    
    # Create plot
    p = surface(xs, ys, zs, 
                color=:lightblue, alpha=0.2,
                legend=false, title=title,
                aspect_ratio=:equal,
                camera=(30, 30),
                xlabel="X", ylabel="Y", zlabel="Z",
                xlims=(-1.1,1.1), ylims=(-1.1,1.1), zlims=(-1.1,1.1))
    
    # Add axes
    plot!([-1.1, 1.1], [0, 0], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [-1.1, 1.1], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [0, 0], [-1.1, 1.1], color=:black, linewidth=1, linestyle=:dash)
    
    # Add vector
    plot!([0, x], [0, y], [0, z], 
          color=:red, linewidth=3, arrow=true)
    
    # Add state label
    annotate!(x*1.05, y*1.05, z*1.05, 
              text("|ψ⟩ = cos(θ/2)|0⟩ + e^iϕ sin(θ/2)|1⟩", 8, :red))
    
    # Add basis state labels
    annotate!(1.1, 0, 0, text("|+⟩", 10))
    annotate!(-1.1, 0, 0, text("|-⟩", 10))
    annotate!(0, 1.1, 0, text("|i+⟩", 10))
    annotate!(0, -1.1, 0, text("|i-⟩", 10))
    annotate!(0, 0, 1.1, text("|0⟩", 10))
    annotate!(0, 0, -1.1, text("|1⟩", 10))
    
    return p
end

# Example usage
θ = π/3  # Polar angle (from z-axis)
ϕ = π/4  # Azimuthal angle (in x-y plane)
plot_bloch_state(θ, ϕ, title="Qubit State on Bloch Sphere")