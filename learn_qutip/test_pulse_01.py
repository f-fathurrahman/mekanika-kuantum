import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("dark_background")


import qutip
from qutip import sigmax, sigmaz, basis, mesolve

ω = 10.0
μ01 = 1.0
H0 = -ω/2 * sigmaz()
V = -μ01 * sigmax()


def my_pulse(t, E0=1.0, ωl=10.0, t0=25.0, τ=1.5, ϕ=0.0):
    return E0 * np.cos(ωl * (t - t0) + ϕ) * np.exp(-(t - t0)**2 / (2*τ**2))


t0 = 25.0
ωl = ω
E0 = 1.2
ϕ = 0
τ = 1.5
args_pulse = {"ωl": ωl, "E0": E0}
H = H0 + V * qutip.coefficient(my_pulse, args=args_pulse)

tlist = np.linspace(0.0, 50.0, 1000)
Epulse = my_pulse(tlist, **args_pulse)
psi0 = basis(2, 0)
P0 = basis(2, 0).proj()
P1 = basis(2, 1).proj()
result = mesolve(H, psi0, tlist, e_ops=[P0, P1], options={"store_states": True})

plt.clf()
plt.plot(tlist, Epulse)
plt.show()

plt.clf()
plt.plot(tlist, result.expect[0])
plt.plot(tlist, result.expect[1])
