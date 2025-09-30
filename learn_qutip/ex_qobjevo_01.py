# Basic setup
import numpy as np
import qutip

Nsize = 4
t = 1.0

# Some basic operators
a = qutip.destroy(Nsize)
ad = qutip.create(Nsize)
n = qutip.num(Nsize)

# Identity
Id = qutip.qeye(Nsize)

# Three types of scalar time dependence: function, string, array
#
# A(t) = ∑_k f_k(t) A_k
#
# f_k(t) are time-dependet scalars
# A_k are constant Qobj objects
#
# Can be represented as list:
# [ A0, [A1, f1], [A2, f2], ... ]
#
# There are various ways to represent fk

# Alternatively, QobjEvo also can be created by multiplication of Qobj with wrapped time dependencies

# Example of constant form, using number operator
constant_form = qutip.QobjEvo([n])

# Function for time-dependence
def my_cos_t(t):
    return np.cos(t)

# We use Qobj multiplied by qutip.coefficient, Qobj will be promoted to QobjEvo
function_form = n + (a + ad)*qutip.coefficient(my_cos_t)

# for something with memory (?)
# SKIPPED


# Evaluation
print("constant_form(2.0) = ", constant_form(2.0))
print("function_form(2.0) = ", function_form(2.0))

