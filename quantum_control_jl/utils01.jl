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

