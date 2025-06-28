import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Kraus, SuperOp
from qiskit.visualization import plot_histogram
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Statevector,DensityMatrix,state_fidelity,partial_trace, Operator
from matplotlib import pyplot as plt
from functools import reduce
from scipy.linalg import expm
import pandas as pd

# Import from Qiskit Aer noise module
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    ReadoutError,
    depolarizing_error,
    pauli_error,
    thermal_relaxation_error,
)


def construct_hamiltonian(N, J, h):
    # definitions
    sigma_x = np.array([[0, 1], [1, 0]])
    sigma_z = np.array([[1, 0], [0, -1]])
    identity = np.eye(2)
    #
    H = np.zeros((2**N, 2**N))
    # Interaction term: -J * sum(sigma_z[i] * sigma_z[i+1])
    for i in range(N):
        for j in range(i+1, N):
            term = 1
            for k in range(N):
                if k == i:
                    term = np.kron(term, sigma_z)
                elif k == j:
                    term = np.kron(term, sigma_z)
                else:
                    term = np.kron(term, identity)
            H += -J * term
    #
    # Transverse field term: -h * sum(sigma_x[i])
    for i in range(N):
        term = 1
        for j in range(N):
            if j == i:
                term = np.kron(term, sigma_x)
            else:
                term = np.kron(term, identity)
        H += -h * term
    #
    return H

def getExactState(n_qubits, J, h, time):
    
    print("n_qubits = ", n_qubits)

    hamiltonian = construct_hamiltonian(n_qubits, J, h)

    # Exact time evolution using matrix exponentiation
    exact_operator = expm(-1j * time * hamiltonian)
    
    # Create an initial state (|00...0>)
    initial_state = np.zeros(2**n_qubits)
    initial_state[0] = 1
    
    # Apply exact time evolution
    pure_state = exact_operator @ initial_state
    pure_state = DensityMatrix(pure_state)

    return pure_state


list_of_epsilon = [i * 0.005 for i in range(21)]
t = 1
J = 1
h = 1
pure_state = getExactState(2, J, h, t)

# Nmax=5
#N = 2
#list_of_purified_fidelity_flag0_n1_noisy=[[] for _ in range(N)]
#list_of_purified_fidelity_flag0_n2_noisy=[[] for _ in range(N)]
#list_of_fidelity=[[] for _ in range(N)]

#for j in range(N): # loop over no. of qubits
#    for i in list_of_epsilon:
#        list_of_purified_fidelity_flag0_n1_noisy[j].append(state_fidelity(getPurifiedRhoWithNoisySWAP1213GHZ(i,t, J, h, j+1,1,0,i,1), pure_state))
#        list_of_purified_fidelity_flag0_n2_noisy[j].append(state_fidelity(getPurifiedRhoWithNoisySWAP1213GHZ(i,t, J, h, j+1,2,0,i,1), pure_state))
#        list_of_fidelity[j].append(state_fidelity(getInputRho(i,t,j+1,1), pure_state))