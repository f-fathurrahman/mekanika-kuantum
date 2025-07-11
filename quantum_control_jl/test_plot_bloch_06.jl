using Plots
using LinearAlgebra
gr()  # Using GR backend for fast static plots

function plot_bloch_trajectory(points::Vector{Vector{Float64}}; 
                              title="Qubit Dynamics on Bloch Sphere",
                              elevation=30, azimuth=30,
                              point_size=4,
                              trajectory_color=:red,
                              initial_color=:green,
                              final_color=:blue)
    
    # Normalize all points to ensure they're on the sphere
    normalized_points = [p / norm(p) for p in points]
    
    # Extract coordinates for plotting
    x = [p[1] for p in normalized_points]
    y = [p[2] for p in normalized_points]
    z = [p[3] for p in normalized_points]
    
    # Create sphere surface
    u = range(0, 2π, length=50)
    v = range(0, π, length=50)
    x_sphere = [cos(ϕ)*sin(θ) for θ in v, ϕ in u]
    y_sphere = [sin(ϕ)*sin(θ) for θ in v, ϕ in u]
    z_sphere = [cos(θ) for θ in v, ϕ in u]
    
    # Create plot with specified view angle
    p = surface(x_sphere, y_sphere, z_sphere, 
               color=:lightblue, alpha=0.2,
               legend=false, title=title,
               aspect_ratio=:equal,
               camera=(elevation, azimuth),
               xlabel="X", ylabel="Y", zlabel="Z",
               xlims=(-1.1,1.1), ylims=(-1.1,1.1), zlims=(-1.1,1.1))
    
    # Add axes
    plot!([-1.1, 1.1], [0, 0], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [-1.1, 1.1], [0, 0], color=:black, linewidth=1, linestyle=:dash)
    plot!([0, 0], [0, 0], [-1.1, 1.1], color=:black, linewidth=1, linestyle=:dash)
    
    # Plot trajectory points
    scatter!(x, y, z, 
            markersize=point_size,
            color=trajectory_color,
            marker=:circle,
            label="Trajectory")
    
    # Highlight initial and final points
    scatter!([x[1]], [y[1]], [z[1]], 
            markersize=point_size+2,
            color=initial_color,
            marker=:star5,
            label="Initial state")
    
    scatter!([x[end]], [y[end]], [z[end]], 
            markersize=point_size+2,
            color=final_color,
            marker=:diamond,
            label="Final state")
    
    # Add basis state labels
    annotate!(1.1, 0, 0, text("|x⟩", 10))
    annotate!(0, 1.1, 0, text("|y⟩", 10))
    annotate!(0, 0, 1.1, text("|z⟩", 10))
    
    return p
end

# Simulate Rabi oscillations (rotation around x-axis)
function rabi_oscillation(times, Ω)
    [begin
        θ = Ω*t
        [0, sin(θ), cos(θ)]
    end for t in times]
end

times = range(0, 2π, length=50)
trajectory = rabi_oscillation(times, 1.0)

# Plot with different view angles
p1 = plot_bloch_trajectory(trajectory, title="Rabi Oscillation - Standard View")
p2 = plot_bloch_trajectory(trajectory, title="Rabi Oscillation - Top View", elevation=90)
p3 = plot_bloch_trajectory(trajectory, title="Rabi Oscillation - Side View", elevation=0)

plot(p1, p2, p3, layout=(1,3), size=(1200,400))


#=
# Simulate evolution under Hadamard gate (rotation around (x+z) axis)
function hadamard_evolution(times)
    [begin
        θ = π*t/2  # Complete Hadamard at t=1
        [sin(θ)/√2, 0, cos(θ)/√2 + 1/√2]
    end for t in times]
end

times = range(0, 1, length=30)
trajectory = hadamard_evolution(times)

plot_bloch_trajectory(trajectory, 
                     title="Hadamard Gate Evolution",
                     trajectory_color=:purple,
                     point_size=5)



function plot_bloch_trajectory_with_lines(points; kwargs...)
    p = plot_bloch_trajectory(points; kwargs...)
    
    # Extract coordinates
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    z = [p[3] for p in points]
    
    # Add connecting lines
    plot!(x, y, z, 
         color=kwargs[:trajectory_color] ? kwargs[:trajectory_color] : :red,
         linewidth=1, alpha=0.5,
         label=false)
    
    return p
end

# Example usage with lines
times = range(0, π, length=20)
trajectory = [[sin(t), 0, cos(t)] for t in times]
plot_bloch_trajectory_with_lines(trajectory, 
                               title="Qubit Trajectory with Connecting Lines",
                               trajectory_color=:orange)
=#