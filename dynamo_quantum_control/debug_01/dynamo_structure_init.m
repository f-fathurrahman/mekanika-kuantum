function dynamo_structure_init(num_spins)

global OC;

disp('<div> ENTER dynamo_structure_init');

OC = struct();

OC.config.uInitial = NaN(2^num_spins);
OC.config.uFinal = NaN(2^num_spins);

OC.config.numControls = NaN;

OC.config.totalTime = NaN;

OC.config.controlsInitialScaling = NaN;

OC.config.normType = 'PSU'; % 'SU' and 'PSU' supported

OC.timeSlots = struct();
OC.timeSlots.nTimeSlots = NaN;
OC.timeSlots.expmFunc = @expm;
OC.timeSlots.calcPfromHfunc = @calcPfromH_expm;

OC.config.gradientFunc = @gradient_exact;
OC.timeSlots.calcPfromHfunc = @calcPfromH_exact_gradient;

OC.const = struct();
OC.const.termination_reason = struct( ...
    'goal_achieved',    1, ...
    'loop_count',       2, ...
    'wall_time',        3, ...
    'cpu_time',         4, ...
    'gradient_norm',    5);
OC.const.termination_reason_str = { ...
    'Goal achieved', ...
    'Loop count limit reached', ...
    'Wall time limit reached', ...
    'CPU time limit reached', ...
    'Minimal gradient norm reached'};

OC.config.hamDrift = NaN(2^num_spins);
OC.config.hamControl = {NaN(2^num_spins), NaN(2^num_spins)};
OC.config.normFunc = @PSU_norm;
OC.config.gradientNormFunc = @gradient_PSU_norm;

disp(OC.config.hamDrift);
disp(OC.config.hamControl);

disp('</div> EXIT dynamo_structure_init');
