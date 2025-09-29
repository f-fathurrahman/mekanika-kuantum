import matplotlib.pyplot as plt
import numpy as np
import qutip
from qutip import Bloch, QobjEvo, basis, sesolve, sigmay, sigmaz

psi = (2.0*basis(2, 0) + basis(2,1)).unit()
b = Bloch()
b.add_states(psi)
b.show()
plt.show()

# Simulation with constant magnetic field
H = sigmaz()
tgrid = np.linspace(0, 10.0, 100)
res = sesolve(H, psi, tgrid, [sigmay()])

plt.plot(tgrid, res.expect[0])
plt.xlabel("Time"), plt.ylabel("<sigma_y>")
plt.show()

# Simulation without computing expectations
res = sesolve(H, psi, tgrid)
b = Bloch()
b.add_states(res.states[1:30])
b.show()

# Simulation with varying magnetic field
def linear(t, args):
    return 0.3 * t

def periodic(t, args):
    return np.cos(0.5 * t)

# Define QobjEvos
H_lin = QobjEvo([[sigmaz(), linear]], tlist=tgrid)
H_per = QobjEvo([[sigmaz(), periodic]], tlist=tgrid)

result_lin = sesolve(H_lin, psi, tgrid, e_ops=[sigmay()])
result_per = sesolve(H_per, psi, tgrid, e_ops=[sigmay()])

# Plot <sigma_y> for linear increasing field strength
plt.plot(tgrid, result_lin.expect[0])
plt.xlabel("Time"), plt.ylabel("<sigma_y>")
plt.show()

plt.plot(times, result_per.expect[0])
plt.xlabel("Time"), plt.ylabel("<sigma_y>")
plt.show()