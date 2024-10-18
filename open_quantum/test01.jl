using MyOpenQuantumTools
using OrdinaryDiffEq

# Define Hamiltonian
# s is dimensionless time s=s/tf where tf is final simulation time
H = DenseHamiltonian(
    [(s)->1-s, (s)->s],
    [σx, σz],
    unit=:ħ
)
# The Hamiltonian is sum of two terms
# H = (1-s)*σx + s*σz

# Final time of simulation
tf = 20

# Initial condition
u0 = PauliVec[1][2]

annealing = Annealing(H, u0)
sol = solve_schrodinger(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6)
