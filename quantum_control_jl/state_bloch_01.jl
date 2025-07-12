using LinearAlgebra

"""
    state_to_bloch(ψ::Vector{ComplexF64}) -> Vector{Float64}

Convert a qubit state vector to Bloch sphere coordinates.
The state vector should be normalized (∑|ψᵢ|² = 1).
"""
function state_to_bloch(ψ::Vector{ComplexF64})
    # Check input is a valid qubit state
    if length(ψ) != 2
        error("State vector must have length 2 for qubit")
    end
    
    # Verify normalization
    if !isapprox(norm(ψ), 1.0, atol=1e-8)
        @warn "State vector not normalized (norm = $(norm(ψ))), normalizing..."
        ψ = ψ / norm(ψ)
    end
    
    # Extract components
    α, β = ψ
    
    # Calculate Bloch coordinates
    x = 2 * real(α * conj(β))
    y = 2 * imag(α * conj(β))
    z = abs2(α) - abs2(β)
    
    return [x, y, z]
end

"""
    bloch_to_state(r::Vector{Float64}) -> Vector{ComplexF64}

Convert Bloch sphere coordinates to a qubit state vector (returns normalized state).
"""
function bloch_to_state(r::Vector{Float64})
    # Check input is valid Bloch vector
    if length(r) != 3
        error("Bloch vector must have length 3")
    end
    
    # Verify normalization
    norm_r = norm(r)
    if !isapprox(norm_r, 1.0, atol=1e-8)
        @warn "Bloch vector not normalized (norm = $norm_r), normalizing..."
        r = r / norm_r
    end
    
    θ = acos(r[3])  # Polar angle
    ϕ = atan(r[2], r[1])  # Azimuthal angle
    
    # Return state vector in standard form
    return [cos(θ/2), exp(im*ϕ) * sin(θ/2)]
end

# Convert |+⟩ state (x=1)
ψ_plus = [1/sqrt(2), 1/sqrt(2)]  # |+⟩ = (|0⟩ + |1⟩)/√2
r_plus = state_to_bloch(ψ_plus)  # Returns [1.0, 0.0, 0.0]

# Convert |-⟩ state (x=-1)
ψ_minus = [1/sqrt(2), -1/sqrt(2)]  # |-⟩ = (|0⟩ - |1⟩)/√2
r_minus = state_to_bloch(ψ_minus)  # Returns [-1.0, 0.0, 0.0]

# Convert |i+⟩ state (y=1)
ψ_iplus = [1/sqrt(2), im/sqrt(2)]  # |i+⟩ = (|0⟩ + i|1⟩)/√2
r_iplus = state_to_bloch(ψ_iplus)  # Returns [0.0, 1.0, 0.0]

# Convert arbitrary state
ψ = normalize([1, 2im])  # Normalized complex state
r = state_to_bloch(ψ)    # Returns Bloch vector

# Simulate unitary evolution of a qubit state
function simulate_unitary(ψ₀::Vector{ComplexF64}, U::Matrix{ComplexF64}, steps::Int)
    [state_to_bloch(U^n * ψ₀) for n in 0:steps-1]
end

# Example: Rabi oscillation between |0⟩ and |1⟩
ψ₀ = [1.0, 0.0]  # |0⟩ state
θ = π/8  # Rotation angle per step
U = [cos(θ) -im*sin(θ); -im*sin(θ) cos(θ)]  # Rotation around x-axis

# Generate trajectory
trajectory = simulate_unitary(ψ₀, U, 16)

# Plot with Qutip style
qutip_style_bloch(trajectory, 
                 title="Unitary Evolution Trajectory",
                 trajectory_color=:red)

# |0⟩ state (north pole)
ψ0 = [1.0, 0.0]
r0 = state_to_bloch(ψ0)  # [0.0, 0.0, 1.0]

# |1⟩ state (south pole)
ψ1 = [0.0, 1.0]
r1 = state_to_bloch(ψ1)  # [0.0, 0.0, -1.0]

# Superposition with phase
ψ = normalize([1.0, exp(im*π/4)])  # Arbitrary state
r = state_to_bloch(ψ)