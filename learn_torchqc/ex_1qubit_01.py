import numpy as np
import torch

from torchqc.states import QuantumState
from torchqc.operators import Operator
from torchqc.dynamics import TDSE

# Setup initial condition
n = 2
basis_states = QuantumState.basis(n)
state = basis_states[0] + 2.0*basis_states[1] # |0> + 2|1>
state.normalize()
# try to play with initial state

# Simulate dynamics with TDSE
T = 10
Δt = 0.01
#time = np.arange(0, T + Δt, Δt, dtype = np.float64) # XXX why float32?
time = np.arange(0, T + Δt, Δt)
matrices = torch.from_numpy(
    np.array([
        np.array([[0, 1], [1, 0]], dtype=np.complex128) for _ in time
    ])
)
hamiltonian = Operator(2, matrices)

states = TDSE(state, hamiltonian, time, Δt)
# states is a list of QuantumState

import matplotlib.pyplot as plt
populations = np.array([
    (torch.abs(state.state_tensor)**2).numpy() for state in states
])

fig, ax = plt.subplots()
ax.plot(time, populations[:,0], label = "P1")
ax.plot(time, populations[:,1], label = "P2")
ax.legend()
plt.grid(True)
plt.savefig("IMG_1qubit_01.png", dpi=150)
plt.show()
