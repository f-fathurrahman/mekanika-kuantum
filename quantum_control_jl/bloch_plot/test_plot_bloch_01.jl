using Plots
using LinearAlgebra

function plot_bloch_vector(r::Vector{Float64}; title="")
    # Ensure the vector is normalized
    norm_r = norm(r)
    if !isapprox(norm_r, 1.0, atol=1e-3)
        @warn "Vector is not normalized (norm = $norm_r), normalizing..."
        r = r / norm_r
    end
    
    # Create sphere surface
    u = range(0, 2π, length=50)
    v = range(0, π, length=50)
    x = [cos(ϕ)*sin(θ) for θ in v, ϕ in u]
    y = [sin(ϕ)*sin(θ) for θ in v, ϕ in u]
    z = [cos(θ) for θ in v, ϕ in u]
    
    # Create plot
    p = surface(x, y, z, 
                color=:lightblue, alpha=0.3, 
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
    plot!([0, r[1]], [0, r[2]], [0, r[3]], 
          color=:red, linewidth=3, arrow=true)
    
    # Add labels
    annotate!(1.1, 0, 0, text("|x⟩", 10))
    annotate!(0, 1.1, 0, text("|y⟩", 10))
    annotate!(0, 0, 1.1, text("|z⟩", 10))
    
    return p
end

# Example usage
r = [1/√3, 1/√3, 1/√3]  # Equal superposition
plot_bloch_vector(r, title="Bloch Sphere Representation")