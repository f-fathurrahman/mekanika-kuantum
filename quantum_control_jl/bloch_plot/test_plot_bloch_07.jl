using Plots
using LinearAlgebra
gr()  # Using GR backend

function qutip_style_bloch(points::Vector{Vector{Float64}};
                          title="",
                          elevation=30, azimuth=30,
                          point_size=4,
                          trajectory_color=:red,
                          initial_color=:green,
                          final_color=:blue,
                          show_axes=true,
                          show_grid=true)

    # Normalize all points
    normalized_points = [p / norm(p) for p in points]
    x, y, z = ([p[i] for p in normalized_points] for i in 1:3)

    # Create sphere with Qutip-style appearance
    u = range(0, 2π, length=100)
    v = range(0, π, length=50)
    
    # Main sphere (semi-transparent with cool/warm color gradient)
    p = surface([cos(ϕ)*sin(θ) for θ in v, ϕ in u],
                [sin(ϕ)*sin(θ) for θ in v, ϕ in u],
                [cos(θ) for θ in v, ϕ in u],
                color=:deep, # ???
                #alpha=0.25,  # Semi-transparent
                legend=false,
                title=title,
                aspect_ratio=:equal,
                camera=(elevation, azimuth),
                xlims=(-1.1,1.1), ylims=(-1.1,1.1), zlims=(-1.1,1.1),
                grid=false,      # We'll add our own grid lines
                background_color=:white)  # White background like Qutip

    # Add Qutip-style grid lines (more prominent than default)
    if show_grid
        for ϕ in range(0, 2π, length=8)[1:end-1]  # Longitudinal lines
            plot!([cos(ϕ)*sin(θ) for θ in range(0, π, length=30)],
                  [sin(ϕ)*sin(θ) for θ in range(0, π, length=30)],
                  [cos(θ) for θ in range(0, π, length=30)],
                  color=:black, linewidth=0.5, alpha=0.3, label=false)
        end
        
        for θ in range(0, π, length=6)[2:end-1]  # Latitudinal lines
            plot!([cos(ϕ)*sin(θ) for ϕ in range(0, 2π, length=50)],
                  [sin(ϕ)*sin(θ) for ϕ in range(0, 2π, length=50)],
                  fill(cos(θ), 50),
                  color=:black, linewidth=0.5, alpha=0.3, label=false)
        end
    end

    # Add Qutip-style axes
    if show_axes
        axis_col = RGB(0.5,0.5,0.5)  # Gray color for axes
        plot!([-1.1, 1.1], [0, 0], [0, 0], color=axis_col, linewidth=1.5, label=false)
        plot!([0, 0], [-1.1, 1.1], [0, 0], color=axis_col, linewidth=1.5, label=false)
        plot!([0, 0], [0, 0], [-1.1, 1.1], color=axis_col, linewidth=1.5, label=false)
        
        # Add axis labels in Qutip style
        annotate!(1.15, 0, 0, text("x", 10, :right, axis_col))
        annotate!(0, 1.15, 0, text("y", 10, :top, axis_col))
        annotate!(0, 0, 1.15, text("z", 10, :top, axis_col))
        annotate!(0, 0, -1.15, text("-z", 10, :bottom, axis_col))
    end

    # Plot trajectory with Qutip-style markers
    scatter!(x, y, z,
            markersize=point_size,
            color=trajectory_color,
            marker=:circle,
            markerstrokewidth=0.5,
            label="Trajectory")

    # Initial and final points with Qutip-style markers
    scatter!([x[1]], [y[1]], [z[1]],
            markersize=point_size+2,
            color=initial_color,
            marker=:pentagon,
            markerstrokewidth=1,
            label="Initial")

    scatter!([x[end]], [y[end]], [z[end]],
            markersize=point_size+2,
            color=final_color,
            marker=:star,
            markerstrokewidth=1,
            label="Final")

    # Add connecting lines with transparency
    plot!(x, y, z,
         color=trajectory_color,
         linewidth=1.5,
         alpha=0.6,
         label=false)

    return p
end

# Simulate a quantum rotation
times = range(0, 2π, length=50)
trajectory = [[sin(t)*cos(2t), sin(t)*sin(2t), cos(t)] for t in times]

# Create Qutip-style plot
qutip_style_bloch(trajectory,
                 title="Qutip-style Bloch Sphere",
                 elevation=45, azimuth=30,
                 trajectory_color=RGB(0.9,0.2,0.2),  # Qutip-like red
                 initial_color=RGB(0.2,0.7,0.2),     # Qutip-like green
                 final_color=RGB(0.1,0.3,0.8))       # Qutip-like blue

#=
function add_qutip_style_features!(plot)
    # Remove axes and ticks
    plot.attr[:xaxis][:showaxis] = false
    plot.attr[:yaxis][:showaxis] = false
    plot.attr[:zaxis][:showaxis] = false
    
    # Add the "floor shadow" like Qutip
    ϕ = range(0, 2π, length=100)
    plot!([cos(θ) for θ in ϕ], [sin(θ) for θ in ϕ], fill(-1.1, 100),
         color=:black, linewidth=1, alpha=0.1, label=false)
    
    # Adjust lighting to match Qutip's shading
    plot.attr[:series_list][1][:extra_kwargs][:light_position] = (1,1,1)
end

# Usage:
p = qutip_style_bloch(trajectory)
add_qutip_style_features!(p)
display(p)
=#