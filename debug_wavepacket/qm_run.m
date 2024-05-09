%% Calculate bound states of the Morse oscillator
qm_setup(); qm_init(); qm_bound(); qm_cleanup();

%% Calculate FC spectrum from overlaps of ground state with all other states
qm_setup(); qm_spec(); aux.qm_fcspec(); qm_cleanup();