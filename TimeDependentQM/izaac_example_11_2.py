# Solve the Schroedinger equation for a wave packet with no initial momentum.

# Solution: Let's use the Runge–Kutta fourth order method for this example,
# using Δx = 0.1, and choosing Δt = 0.01 ≤ Δx^2 for stability of our finite-difference
# discretization.

import numpy as np

# create the x grid
Δx = 0.1
x = np.arange(-10, 10+Δx, Δx)
N = len(x)

# create the initial wave packet
psi0 = np.exp(-x**2/4)/((2*np.pi)**(1/4))

# the potential is zero in this case
V = np.zeros([N])

# construct the 4th order FD matrix
g = -5j/(4*Δx**2) - 1j*V
a = 1j/(24*Δx**2)
diag = np.diag(g)
off_diag1 = np.diag([16*a]*(N-1), 1) + np.diag([16*a]*(N-1), -1)
off_diag2 = np.diag([-a]*(N-2), 2) + np.diag([-a]*(N-2), -2)
M = diag + off_diag1 + off_diag2

# create the time grid
dt = 0.01
t = np.arange(0, 20+dt, dt)
steps = len(t)
print("steps = ", steps)
# create an array containing wavefunctions for each step
y = np.zeros([steps, N], dtype=np.complex128)
y[0] = psi0

# the RK4 method
for i in range(0, steps-1):
    k1 = np.dot(M, y[i])
    k2 = np.dot(M, y[i] + k1*dt/2)
    k3 = np.dot(M, y[i] + k2*dt/2)
    k4 = np.dot(M, y[i] + k3*dt)
    y[i+1] = y[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6


import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("dark_background")
matplotlib.rcParams.update({
    "axes.grid": True,
    "grid.color": "grey"
})

for idx_time in [0, 100, 200, 1000, 1500]:
    plt.plot(x, np.real(y[idx_time]))
plt.title("Real part")
plt.show()

for idx_time in [0, 100, 200, 1000, 1500]:
    plt.plot(x, np.imag(y[idx_time]))
plt.title("Imag part")
plt.show()

for idx_time in [0, 100, 200, 1000, 1500]:
    plt.plot(x, np.abs(y[idx_time]))
plt.title("Abs part")
plt.show()