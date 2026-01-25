import numpy as np
from scipy import linalg

from qutip import Qobj, Bloch, sigmax, sigmay, sigmaz, expect

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("dark_background")
matplotlib.rcParams.update({
    "axes.grid" : True,
    "grid.color": "gray"
})

def plot_traj_points_bloch(traj_psi):
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(projection="3d"))
    ax.axis("square") # to get a nice circular plot
    Nvecs = len(traj_psi)
    x = np.zeros(Nvecs)
    y = np.zeros(Nvecs)
    z = np.zeros(Nvecs)
    for i in range(Nvecs):
        psi = Qobj(traj_psi[i])
        x[i] = expect(sigmax(), psi)
        y[i] = expect(sigmay(), psi)
        z[i] = expect(sigmaz(), psi)
    #
    for i in range(Nvecs):
        b = Bloch(fig=fig, axes=ax)
        b.add_points([x[:i+1], y[:i+1], z[:i+1]])
        b.add_vectors([x[i], y[i], z[i]])
        b.render()
        plt.savefig(f"IMG_traj_{i:04d}.png", dpi=150)
        print(f"Fig. {i} is done")
    return


def plot_traj_bloch(traj_psi, filename):
    plt.close()
    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(projection="3d"))
    ax.axis("square") # to get a nice circular plot
    Nvecs = len(traj_psi)
    x = np.zeros(Nvecs)
    y = np.zeros(Nvecs)
    z = np.zeros(Nvecs)
    for i in range(Nvecs):
        psi = Qobj(traj_psi[i])
        x[i] = expect(sigmax(), psi)
        y[i] = expect(sigmay(), psi)
        z[i] = expect(sigmaz(), psi)
    #
    b = Bloch(fig=fig, axes=ax)
    b.add_points([x, y, z])
    b.render()
    plt.savefig(filename, dpi=150)
    return


ket0 = np.array([1.0, 0.0])
ket1 = np.array([0.0, 1.0])

Sx = 1/2 * np.array([[0, 1],[ 1, 0]], dtype=np.complex128)
Sy = 1/2 * np.array([[0, -1j],[1j, 0]], dtype=np.complex128)
Sz = 1/2 * np.array([[1, 0],[0, -1]], dtype=np.complex128)

# Return Hamiltonian at time index j
def eval_hamiltonian(Jt):
    H_drift = Sx
    H_control = 4*Jt*Sz
    return H_drift + H_control

T = 2*np.pi       
Ntimes = 200
dt = T/Ntimes

NitersMax = 30
fidelity = []

# This is a projector to |1>
P1 = np.outer(ket1, ket1)

traj_fw = np.zeros((Ntimes+1,2), dtype=np.complex128) # forward trajectory
traj_fw[0] = ket0

seq = np.random.rand(Ntimes) # control
seq_f = np.zeros(Ntimes)
seq_f[:] = seq[:] # copy

# Propagate
for i in range(Ntimes):
    traj_fw[i+1] = linalg.expm(-(1j) * eval_hamiltonian(seq[i]) * dt).dot(traj_fw[i])

traj_bw = np.zeros((Ntimes+1,2), dtype=np.complex128) # backward trajectory
fidelity.append(np.abs(traj_fw[-1].dot(ket1))**2)

traj_bw[-1] = P1 @ traj_fw[-1]
dj = 0.01  # variation in control sequence? for computing dH

for i in range(NitersMax):
    for j in reversed(range(Ntimes)):
        traj_bw[j] = linalg.expm((1j) * eval_hamiltonian(seq[j]) * dt) @ traj_bw[j+1]
    for k in range(Ntimes):
        dH = (eval_hamiltonian(seq[k]+dj) - eval_hamiltonian(seq[k]-dj)) / (2*dj)
        # Update sequence for time k
        dseq = traj_bw[k].conj().T @ dH @ traj_fw[k]
        #print("dseq = ", dseq)
        seq_f[k] = seq[k] + np.imag(traj_bw[k].conj().T @ dH @ traj_fw[k])
        # forward trajectory
        traj_fw[k+1] = linalg.expm(-(1j) * eval_hamiltonian(seq_f[k]) * dt) @ traj_fw[k]
        seq[:] = seq_f[:]
    plot_traj_bloch(traj_fw, f"IMG_traj_iter_{i+1:04d}.png")
    fidelity.append(np.abs(traj_fw[-1].dot(ket1))**2)
    #fidelity[i+1] += (np.absolute(psi[-1,-1]))**2
    print(f"i={i:4d}, fidelity = {fidelity[i+1]:.10f}")
    if 1 - fidelity[i+1] < 1e-8:
        print("Converged")
        break
    traj_bw[-1] = P1 @ traj_fw[-1]

print('final_fidelity = ', fidelity[-1])
plt.clf()
plt.plot(seq)
plt.savefig("IMG_control_01.png", dpi=150)

