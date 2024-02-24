clear all; close all;

% Global variables are declared here
global hamilt space state;

disp('***********************************')
disp('DEBUG Morse oscillator (OH radical)')
disp('***********************************')

% Initialize state (using wavefunction)
state = wave(); % will call a constructor
state.save_export = false;

% Setup Hamiltonian
hamilt.coupling = ham.coupling;
hamilt.eigen    = ham.eigen;
hamilt.truncate = ham.truncate;

% Spatial discretization
space.dof{1} = dof.fft;  % using fft grid
space.dof{1}.mass = 1728.539; % Reduced mass
space.dof{1}.n_pts = 128;  % Number of grid points
space.dof{1}.x_min = 1.0; % Lower bound of grid 
space.dof{1}.x_max = 10.0; % Upper bound of grid
space.dof{1}.periodic = false; % without PBC

% Hamiltonian operator option (what these parameters affects?)
hamilt.truncate.e_min  =  0.0; % Lower truncation of energy
hamilt.truncate.e_max  =  1.0; % Upper truncation of energy

% See +pot.morse constructor
hamilt.pot{1,1} = pot.morse;  % Morse potential (this will call pot.morse constructor)
hamilt.pot{1,1}.d_e  = 0.1994; % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821; % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189; % Range parameter

% Select eigen/values/functions
hamilt.eigen.start = 0;
hamilt.eigen.stop  = 2;
