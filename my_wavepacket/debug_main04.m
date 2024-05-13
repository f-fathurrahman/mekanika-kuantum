clear all; close all; clc

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

% Initialize state (using wavefunction)
state = WaveClass(); % will call a constructor (defined in @wave/wave)
state.save_export = false;
% state object is initialized here

% Calling a few constructors: temporal discretization
time_var.steps  = tmp.StepsClass();
time_var.efield = tmp.EfieldClass();

% Setup Hamiltonian
hamilt.coupling = CouplingClass();
hamilt.eigen    = EigenClass();
hamilt.truncate = TruncateClass();

prt.disp('***************************************')
prt.disp('Squeezed state of a harmonic oscillator')
prt.disp('***************************************')

% Spatial discretization
space.dof{1}       = FFTClass();            % Use Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -10;                % Lower bound of the grid
space.dof{1}.x_max = 10;                 % Upper bound of the grid
space.dof{1}.mass  = 1;                  % Mass for the kinetic energy

% Temporal discretization
time_var.steps.m_start  = 0;               % Index of initial time step
time_var.steps.m_stop   = 2;               % Index of final time step
time_var.steps.m_delta  = pi/10;             % Size of time steps 
time_var.steps.s_number = 100;              % Number of sub steps per time step

% Initial wave function
time_var.dof{1}       = init.GaussianClass();          % Gaussian-shaped wavepacket
time_var.dof{1}.width = sqrt(2);             % Width 
time_var.dof{1}.pos_0 = -3.0;                % Center in position representation
time_var.dof{1}.mom_0 =  0.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  = 50.0;           % Upper truncation of energy

hamilt.pot{1,1}        = pot.TaylorClass();     % Taylor expansion
hamilt.pot{1,1}.coeffs = [0;1];          % Force constant

%
% .... Plotting stuffs are removed
%


%
% Steps of qm_propa are started here ....
%

% Set propagator
time_var.propa = tmp.wave.strang;
