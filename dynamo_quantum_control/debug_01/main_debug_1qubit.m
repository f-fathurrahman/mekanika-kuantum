% What this demo does

clear variables; close all;

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

intialize_timeslot_controls(initial_controls);

% Which time slots do you want to modify
controls_mask = true(OC.timeSlots.nTimeSlots, OC.config.numControls);

controls_mask(1:100) = false;
initial_controls(1:100) = 0.0;

%controls_mask = false(OC.timeSlots.nTimeSlots, OC.config.numControls);
%controls_mask(200:300) = true;
%initial_controls(1:199) = 0.0;
%initial_controls(301:500) = 0.0;

% Final preperatory configuration
initSetNorm(); % Calculates the norm of the <initial | final> to scale subsequent norms

% Now do the actual search

termination_conditions = struct( ...
  'loop_count',           1e10, ...
  'goal',                 1 - 1e-6, ...
  'wall_time_to_stop',    180, ...
  'cputime_to_stop',      180, ...
  'gradient_norm_min',    1e-20);

wall0 = now();
cpu0 = cputime();

% Which sort of gradient to use
OC.config.gradientFunc = @gradient_exact;
% When computing exact gradient, we get exponentiation for free due to the
% eigendecomposition (see paper for details)
OC.timeSlots.calcPfromHfunc = @calcPfromH_exact_gradient;

OC.config.BFGS = struct('fminopt', struct('Display', 'off'));

% Here we will do the minimization
% Using GRAPE
%termination_reason = BFGS_search_function(controls_mask, termination_conditions);
% Or Krotov
termination_reason = Krotov_search_function(controls_mask, termination_conditions);

fprintf('Fidelity reached: 1 - %g\n', 1-get_current_value())
fprintf('Wall time: %g\n', (now()-wall0)*(24*60*60))
fprintf('CPU time:  %g\n', cputime()-cpu0)
fprintf('Termination reason: %s\n\n\n', OC.const.termination_reason_str{termination_reason});

plot(cumsum(OC.timeSlots.tau), OC.timeSlots.currPoint.controls);
title('Optimized Control Sequences');
xlabel('Time');
ylabel('Control amplitude');
grid on;
axis tight;
