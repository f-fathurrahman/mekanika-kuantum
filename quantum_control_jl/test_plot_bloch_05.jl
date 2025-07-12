using LinearAlgebra
using Plots
gr()  # Set GR backend explicitly

function bloch_sphere_view(r::Vector{Float64}; 
                          title="Bloch Sphere",
                          elevation=30, azimuth=30)
    # Normalize input vector
    r = r / norm(r)
    
    # Sphere surface
    u = range(0, 2π, length=50)
    v = range(0, π, length=50)
    x = [cos(ϕ)*sin(θ) for θ in v, ϕ in u]
    y = [sin(ϕ)*sin(θ) for θ in v, ϕ in u]
    z = [cos(θ) for θ in v, ϕ in u]
    
    # Create plot with specified view angle
    p = surface(x, y, z, 
               color=:lightblue, alpha=0.3,
               legend=false, title=title,
               aspect_ratio=:equal,
               camera=(elevation, azimuth),  # Control view here
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

# Standard view
r = [1/√3, 1/√3, 1/√3]
p1 = bloch_sphere_view(r, title="Standard View")

# View from top down (looking at z-axis)
p2 = bloch_sphere_view(r, title="Top View", elevation=90, azimuth=0)

# View along x-axis
p3 = bloch_sphere_view(r, title="X-Axis View", elevation=0, azimuth=0)

# View along y-axis
p4 = bloch_sphere_view(r, title="Y-Axis View", elevation=0, azimuth=90)

# Display all views
plot(p1, p2, p3, p4, layout=(2,2), size=(800,800))

#=
# Create rotation animation
anim = @animate for azim in range(0, 360, length=36)
    bloch_sphere_view([1/√2, 1/√2, 0], 
                     title="Rotating Bloch Sphere",
                     azimuth=azim,
                     elevation=30)
end
gif(anim, "bloch_rotation.gif", fps=10)
=#