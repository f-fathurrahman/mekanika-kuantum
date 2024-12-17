using MyQuantumDynamics

dt = 0.125
ntimes = 200

ρ0 = [1.0+0.0im 0.0; 0.0 0.0]

# two levels Hamiltonian
H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)

times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes);

