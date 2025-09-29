import matplotlib.pyplot as plt
import numpy as np
import qutip
from qutip import Qobj, Bloch, QobjEvo, basis, sesolve, sigmay, sigmaz, sigmax


psi = (2.0*basis(2, 0) + basis(2,1)).unit()

el_chg = 1.0
ħ = 1.0
m = 1.0
B0 = 10.0
B1 = 1.0
ωL = B0*el_chg/m
ω1 = ωL

# Simulation with constant magnetic field
H0 = sigmaz()
tgrid = np.linspace(0, 10.0, 100)

def periodic_field(t):
    el_chg = 1.0
    ħ = 1.0
    m = 1.0
    B0 = 10.0
    B1 = 1.0
    ωL = B0*el_chg/m
    ω1 = ωL
    return el_chg*ħ*B1/(2*m) * np.cos(ω1*t)

H_per = -0.5*B0*sigmaz() + el_chg*ħ*B1/(2*m) * sigmax() * qutip.coefficient(periodic_field)

ket0 = basis(2, 0).data_as("ndarray")
ket1 = basis(2, 1).data_as("ndarray")
P0 = qutip.Qobj(np.outer(ket0, ket0))
P1 = qutip.Qobj(np.outer(ket1, ket1))
result_per = sesolve(H_per, psi, tgrid, e_ops=[P0, P1])

"""
# Define QobjEvos
H_per = QobjEvo([ [sigmaz(), periodic_field] ], tlist=tgrid)
result_per = sesolve(H_per, psi, tgrid, e_ops=[sigmay()])
"""