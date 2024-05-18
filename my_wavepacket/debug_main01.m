clear all; close all;

% Let's remove these global variables
% global hamilt space state info;

% Initializes general information and sets up log files.
% Running under Matlab or Octave
if exist ('OCTAVE_VERSION', 'builtin')
    info.system = 'Octave';
    warning("off", "Octave:shadowed-function");
    pkg load statistics
    pkg load control
else
    info.system = 'Matlab';
end

% XXX: This is currently disabled as it uses and modify global variable heavily.
%prt.init(mfilename('fullpath'), info);

prt.disp('***********************************')
prt.disp('DEBUG Morse oscillator (OH radical)')
prt.disp('***********************************')

% Initialize state (using wavefunction)
state = WaveClass(); % will call a constructor (defined in @wave/wave)
state.save_export = false;
% state object is initialized here

% Calling a few constructors: temporal discretization
% Not really needed for bound states calculation
time_var.steps  = tmp.StepsClass;
time_var.efield = tmp.EfieldClass;

% Setup Hamiltonian
hamilt.coupling = CouplingClass();
hamilt.eigen    = EigenClass();
hamilt.truncate = TruncateClass();
% In Matlab, there is no need to initialize hamilt structure


space.n_dim = 1;
% Initialize space.dof{1} as fft object
space.dof{1} = FFTClass();  % using fft grid
space.dof{1}.mass = 1728.539; % Reduced mass
space.dof{1}.n_pts = 128;  % Number of grid points
space.dof{1}.x_min = 1.0; % Lower bound of grid
space.dof{1}.x_max = 10.0; % Upper bound of grid
space.dof{1}.periodic = false; % without PBC

% Hamiltonian operator option (what these parameters affects?)
hamilt.truncate.e_min  =  0.0; % Lower truncation of energy
hamilt.truncate.e_max  =  1.0; % Upper truncation of energy

% Initialize potential
hamilt.pot{1} = pot.MorseClass();
hamilt.pot{1}.d_e  = 0.1994; % Dissociation energy
hamilt.pot{1}.r_e  = 1.821; % Equilibrium length
hamilt.pot{1}.alf  = 1.189; % Range parameter

% Select eigen/values/functions
hamilt.eigen.start = 0;
hamilt.eigen.stop  = 2;

% Initialize spatial discretization for each degree of freedom
space = dof_init(state, space);
% This call will modify or update global variable space

% Initialize Hamiltonian operator
[hamilt, space, time_var] = state.init_ham(hamilt, space, time_var);


% initialize the eigen object
% at this point hamiltonian matrix is not yet calculated,
% but its size and several fields are already defined
time_var = hamilt.eigen.init(hamilt, space, time_var);


% Calculate Hamiltonian matrix elements here
hamilt.eigen.setup(hamilt, space);
% After this call, hamilt.eigen.matrix is available

% symmetrization, this is optional, only implemented for 1d
%symm(hamilt.eigen);

%diag(hamilt.eigen);
hamilt.eigen.diag();