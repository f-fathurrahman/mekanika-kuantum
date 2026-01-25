% What this demo does

clear variables; close all;

file_diary = 'LOG_debug_cost_func.html';
file_id = fopen(file_diary, 'w');
fclose(file_id);
diary(file_diary);

% This short demo will optimize a simple two-qubit QFT gate generation problem
% using the DYNAMO package, with the GRAPE search algorithm

% Define the physics of the problem - part 1
nSpins = 1;

% Preparations

dynamo_structure_init(nSpins); % All definitions are in a global variable called OC

global OC; % and now we can access it too

SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
SI = eye(2)/2;

% Define the physics of the problem - part 2
nSpins = 1;

% Drift Hamiltonian
OC.config.hamDrift = -0.5*SZ;
OC.config.hamControl = {SX};

disp('OC.config.hamDrift = ');
disp(OC.con.hamDrift);

disp('OC.config.hamControl = ');
disp(OC.config.hamControl);

OC.config.uInitial = eye(2^nSpins);
OC.config.uFinal = SZ;

% How much time do we have to drive the system?
% The value specified here has empirically been shown to work well
OC.config.totalTime = 5.0; %6*(nSpins-1)/1;

% Set space for goal function
switch 'PSU'
  case 'PSU'
    % "I don't care about global phase"
    OC.config.normFunc = @PSU_norm;
    OC.config.gradientNormFunc = @gradient_PSU_norm;
  case 'SU'
    % "I care about global phase"
    OC.config.normFunc = @SU_norm;
    OC.config.gradientNormFunc = @gradient_SU_norm;
  otherwise
      error ('Currently only SU and PSU norms are supported');
end

% Time-slot configuration.
% Assumption is of equally-sized timeslots, unless otherwise specified in OC.timeSlots.tau
% Number of time slices to specify the control fields
OC.timeSlots.nTimeSlots = 500;
OC.config.controlsInitialScaling = 1;

% Generate random initial values
%randseed(101);
initial_controls = OC.config.controlsInitialScaling * ( ...
    rand(OC.timeSlots.nTimeSlots, length(OC.config.hamControl)) - 0.5);

initialize_timeslot_controls(initial_controls);

% Which time slots do you want to modify
controls_mask = true(OC.timeSlots.nTimeSlots, OC.config.numControls);

%controls_mask(1:100) = false;
%initial_controls(1:100) = 0.0;

%controls_mask = false(OC.timeSlots.nTimeSlots, OC.config.numControls);
%controls_mask(200:300) = true;
%initial_controls(1:199) = 0.0;
%initial_controls(301:500) = 0.0;

% Final preperatory configuration
initSetNorm(); % Calculates the norm of the <initial | final> to scale subsequent norms


update_timeslot_controls(initial_controls, controls_mask);

%v = -get_current_value();
%v = OC.config.normFunc(get_current_value_Phi0_norm());

%grad = -OC.config.gradientNormFunc(subspace_mask);
%last_grad_norm = sum(sum(grad.*grad));

diary off;
