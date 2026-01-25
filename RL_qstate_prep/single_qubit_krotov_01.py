import numpy as np
from scipy import linalg

sx = 1/2 * np.matrix([[0, 1],[ 1, 0]], dtype=complex)
sy = 1/2 * np.matrix([[0, -1j],[1j, 0]], dtype=complex)
sz = 1/2 * np.matrix([[1, 0],[0, -1]], dtype=complex)

# Return Hamiltonian at time index j
def hamiltonian(j):
    J = 4
    H = (j) * J * sz + sx
    return H


T = 2*np.pi       
Ntimes = 200
dt = T/Ntimes

Niters = 500
fidelity = np.zeros(Niters + 1)

# This is a projector to |1>
observable = np.matrix( np.zeros(shape=(2,2), dtype=complex) )
observable[-1, -1] = 1

psi = np.matrix(np.zeros(shape=(2, Ntimes+1), dtype=complex)) # forward trajectory
psi[0,0] = 1
pseudo = np.matrix(np.zeros(shape=(2, Ntimes+1), dtype=complex)) # backward trajectory

seq = np.random.rand(Ntimes) # control
seq_f = np.zeros(Ntimes)
seq_f[:] = seq[:] # copy

# Propagate
for i in range(Ntimes):
    psi[:,i+1] = linalg.expm(-(1j) * hamiltonian(seq[i]) * dt).dot(psi[:,i])
fidelity[0] = (np.absolute(psi[-1,-1]))**2
pseudo[:,-1] = observable.dot(psi[:,-1])
dj = 0.01  # variation in control sequence? for computing dH

for i in range(Niters):
    for j in reversed(range(Ntimes)):
        pseudo[:,j] = linalg.expm((1j) * hamiltonian(seq[j]) * dt).dot(pseudo[:,j+1])
    for k in range(Ntimes):
        dH = (hamiltonian(seq[k]+dj) - hamiltonian(seq[k]-dj)) / (2*dj)
        # Update sequence for time k
        seq_f[k] = seq[k] + ( pseudo[:,k].conj().T.dot(dH.dot(psi[:,k])) ).imag[0,0]
        # forward trajectory
        psi[:,k+1] = linalg.expm(-(1j) * hamiltonian(seq_f[k]) * dt).dot(psi[:,k])
        seq[:] = seq_f[:]
    fidelity[i+1] += (np.absolute(psi[-1,-1]))**2
    print(f"i={i}, fidelity = {fidelity[i+1]}")
    pseudo[:,-1] = observable.dot(psi[:,-1])

print('final_fidelity = ', fidelity[-1])

