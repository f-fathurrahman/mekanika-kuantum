import numpy as np
import goat_quantumcontrol as Qgoat

sigmax = np.array([
    [0, 1],
    [1, 0]
])
#self.sy = [0 -1i; 1i 0];
sigmaz = np.array([
    [1, 0],
    [0, -1]
])

X_gate = sigmax

#-----System parameters------
# define the drift Hamiltonian
H0 = sigmax

# define the control Hamiltonian
Hdrive = sigmaz
# define the target gate
Utarget = X_gate

#-----Pulse parameters------
# define the number of time intervals
n_ts = 1000
# define the evolution time
evo_time = 3
# define the number of amps
num_of_amps = 2

#-----Optimization parameters--
# define the number of maximal iterations
max_iter = 200

# create an instance of the Pulse class to be used
fourier_pulse = Qgoat.FourierPulseWithEnvelope(n_ts=n_ts,
                                               evo_time=evo_time,
                                               num_of_amps=num_of_amps,
                                               window=None)

# create an instance of the Optimizer class
optimizer = Qgoat.Optimizer(H0=H0, Hdrive=Hdrive,
                            target=Utarget,
                            pulse=fourier_pulse,
                            max_iter=max_iter,
                            printProgress=True)

# run the optimization
optimizer.run_optimization()