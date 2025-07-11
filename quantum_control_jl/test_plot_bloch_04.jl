using PlotlyJS

function plotly_bloch_sphere(θ, ϕ; title="Bloch Sphere")
    # Convert to Cartesian coordinates
    x = sin(θ) * cos(ϕ)
    y = sin(θ) * sin(ϕ)
    z = cos(θ)
    
    # Create sphere
    u = range(0, 2π, length=50)
    v = range(0, π, length=50)
    x_surf = [cos(ϕ)*sin(θ) for θ in v, ϕ in u]
    y_surf = [sin(ϕ)*sin(θ) for θ in v, ϕ in u]
    z_surf = [cos(θ) for θ in v, ϕ in u]
    
    # Create traces
    sphere = surface(x=x_surf, y=y_surf, z=z_surf,
                     colorscale=[[0, "lightblue"], [1, "lightblue"]],
                     opacity=0.6, showscale=false)
    
    axes = [
        scatter3d(x=[-1.1, 1.1], y=[0, 0], z=[0, 0], 
                 mode="lines", line=attr(color="black", width=2, dash="dash")),
        scatter3d(x=[0, 0], y=[-1.1, 1.1], z=[0, 0], 
                 mode="lines", line=attr(color="black", width=2, dash="dash")),
        scatter3d(x=[0, 0], y=[0, 0], z=[-1.1, 1.1], 
                 mode="lines", line=attr(color="black", width=2, dash="dash"))
    ]
    
    vector = scatter3d(x=[0, x], y=[0, y], z=[0, z],
                      mode="lines+markers",
                      line=attr(color="red", width=6),
                      marker=attr(size=4, color="red"))
    
    labels = scatter3d(
        x=[1.15, -1.15, 0, 0, 0, 0],
        y=[0, 0, 1.15, -1.15, 0, 0],
        z=[0, 0, 0, 0, 1.15, -1.15],
        mode="text",
        text=["|+⟩", "|-⟩", "|i+⟩", "|i-⟩", "|0⟩", "|1⟩"],
        textfont=attr(size=14)
    )
    
    layout = Layout(
        title=title,
        scene=attr(
            xaxis=attr(title="X (|+⟩)", range=[-1.2, 1.2]),
            yaxis=attr(title="Y (|i+⟩)", range=[-1.2, 1.2]),
            zaxis=attr(title="Z (|0⟩)", range=[-1.2, 1.2]),
            aspectratio=attr(x=1, y=1, z=1)
        ),
        showlegend=false
    )
    
    Plot([sphere, axes..., vector, labels], layout)
end

# Example usage
plotly_bloch_sphere(π/3, π/4)  # θ=π/3, ϕ=π/4