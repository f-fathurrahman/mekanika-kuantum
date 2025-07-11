using Plots
using LinearAlgebra
plotlyjs()  # Set the interactive backend

function interactive_bloch_sphere(r::Vector{Float64}; title="Bloch Sphere")
    # Normalize the vector if needed
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
    
    # Create plot with PlotlyJS backend
    p = surface(x, y, z, 
                color=:lightblue, opacity=0.6,
                legend=false, title=title,
                aspect_ratio=:equal,
                xlabel="X (|+⟩)", ylabel="Y (|i+⟩)", zlabel="Z (|0⟩)",
                xlims=(-1.1,1.1), ylims=(-1.1,1.1), zlims=(-1.1,1.1))
    
    # Add axes
    plot!([-1.1, 1.1], [0, 0], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [-1.1, 1.1], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [0, 0], [-1.1, 1.1], color=:black, linewidth=1, linestyle=:dash)
    
    # Add Bloch vector
    plot!([0, r[1]], [0, r[2]], [0, r[3]], 
          color=:red, linewidth=4, arrow=true)
    
    # Add basis state labels (positioned slightly outside the sphere)
    scatter!([1.15, -1.15, 0, 0, 0, 0], 
             [0, 0, 1.15, -1.15, 0, 0], 
             [0, 0, 0, 0, 1.15, -1.15],
             marker=:circle, markersize=0,  # invisible markers
             series_annotations=["|+⟩", "|-⟩", "|i+⟩", "|i-⟩", "|0⟩", "|1⟩"])
    
    return p
end

# Example usage
r = [1/√2, 1/√2, 0]  # Equal superposition of |0⟩ and |1⟩ with phase
interactive_bloch_sphere(r)