Pkg.activate("quantumtoolbox", shared=true)

using QuantumToolbox

N = 20 # cufoff of the Hilbert space dimension
ω = 1.0 # frequency of the harmonic oscillator

a = destroy(N) # annihilation operator
H = ω * a' * a

γ = 0.1 # damping rate
ψ₀ = fock(N, 3)
tlist = range(0, 10, 100)

c_ops = [sqrt(γ)*a]
e_ops = [a' * a]

sol = mesolve(H, ψ₀, tlist, c_ops, e_ops=e_ops)