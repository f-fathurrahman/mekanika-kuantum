% Initialize (quantum) Hamiltonian operator for
% use with wave functions represented on grids

% function [hamilt, space, time_var] = init_ham(state, hamilt, space, time_var)

%% Close coupling scheme
init(hamilt.coupling, state);
disp(hamilt.coupling);

%% Optionally truncate all energies below/above thresholds    
if ~isfield(hamilt, 'TruncateClass') % If not already existing
    hamilt.truncate = TruncateClass(); % Construct object
end
init(hamilt.truncate);
disp(hamilt.truncate);

%% Grid representation of kinetic energy (diabatic: scalar)

% Required for some reason inside grid_.../init_kin
time_var.steps.s_delta = 1e-10; % set to dummy value
% ffr: disabled

% Kinetic operators associated with each type of DVR
for k = 1:space.n_dim
    %init_kin(space.dof{k}, space, time_var, hamilt, 1);
    space.dof{k}.init_kin(space, time_var, hamilt, 1)
end

% "Extra" kinetic operators
if isfield(hamilt, 'kin')
    for k = 1:length(hamilt.kin)
        init_kin(hamilt.kin{k}, 1);
    end
end

%% Set up potential energy (diabatic: matrix representation)
if ~isfield(hamilt, 'pot')
    hamilt.pot = cell(hamilt.coupling.n_eqs);
end
for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        if isempty (hamilt.pot{m,n})
            hamilt.pot{m,n} = pot.generic;
        end
        hamilt.pot{m,n}.row = m; % Tell each POT its channels
        hamilt.pot{m,n}.col = n; % Tell each POT its channels
        init( hamilt.pot{m,n}, space );
        grid( hamilt.pot{m,n}, space );
    end
end

%% Set up additional multiplicative operators
if isfield(hamilt, 'amo')
    for p = 1:length(hamilt.amo)
        if isempty (hamilt.amo{p})
            hamilt.amo{p} = amo.generic;
        end
        hamilt.amo{p}.ind = p; % Tell each AMO its index
        init ( hamilt.amo{p} );
        disp ( hamilt.amo{p} ); 
        grid ( hamilt.amo{p} );
    end
else
    prt.disp('***************************************************************')
    prt.disp('No additional multiplicative operators available  ')
    prt.disp('***************************************************************')
    prt.disp('   ')
end

%% Set up system-bath coupling (diabatic: matrix representation)
if isfield(hamilt, 'sbc')
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            if isempty (hamilt.sbc{m,n})
                hamilt.sbc{m,n} = sbc.generic;
            end
            hamilt.sbc{m,n}.row = m; % Tell each SBC its channels
            hamilt.sbc{m,n}.col = n; % Tell each SBC its channels
            init ( hamilt.sbc{m,n} );
            disp ( hamilt.pot{m,n} );
            grid ( hamilt.sbc{m,n} );
        end
    end
else
    prt.disp('***************************************************************')
    prt.disp('No system-bath coupling available                 ')
    prt.disp('***************************************************************')
    prt.disp(' ')
end

%% Set up negative imaginary potential (diabatic: vector)
if isfield(hamilt, 'nip')
    for m = 1:hamilt.coupling.n_eqs
        if isempty (hamilt.nip{m})
            hamilt.nip{m} = nip.generic;
        end
        hamilt.nip{m}.ind = m; % Tell each NIP its channel
        init ( hamilt.nip{m} );
        disp ( hamilt.nip{m} );
        grid ( hamilt.nip{m} );
    end
else
    prt.disp('***************************************************************')
    prt.disp('No negative imaginary potentials available        ')
    prt.disp('***************************************************************')
    prt.disp(' ')
end

%% Set up dipole moments (diabatic: matrix representation)
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            for m = 1:hamilt.coupling.n_eqs
                for n = m:hamilt.coupling.n_eqs
                    if isempty (hamilt.dip{p}{m,n})
                        hamilt.dip{p}{m,n} = dip.generic;
                    end
                    hamilt.dip{p}{m,n}.pol = p; % Tell each DIP its polarization
                    hamilt.dip{p}{m,n}.row = m; % Tell each DIP its channels
                    hamilt.dip{p}{m,n}.col = n; % Tell each DIP its channels
                    init ( hamilt.dip{p}{m,n} );
                    disp ( hamilt.dip{p}{m,n} );
                    grid ( hamilt.dip{p}{m,n} );
                end
            end
        end
    end
else
    prt.disp('***************************************************************')
    prt.disp('No dipole moments available                       ')
    prt.disp('***************************************************************')
    prt.disp('   ')
end

%% Set up polariabilities (diabatic: matrix representation)
if isfield(hamilt, 'pol')
    for p = 1:size(hamilt.pol,1)
        for q = p:size(hamilt.pol,2)
            if ~isempty (hamilt.pol{p,q})
                for m = 1:hamilt.coupling.n_eqs
                    for n = m:hamilt.coupling.n_eqs
                        if isempty (hamilt.pol{p,q}{m,n})
                            hamilt.pol{p,q}{m,n} = pol.generic;
                        end
                        hamilt.pol{p,q}{m,n}.p_1 = p; % Tell each POL its polarization
                        hamilt.pol{p,q}{m,n}.p_2 = q; % Tell each POL its polarization
                        hamilt.pol{p,q}{m,n}.row = m; % Tell each POL its channels
                        hamilt.pol{p,q}{m,n}.col = n; % Tell each POL its channels
                        init ( hamilt.pol{p,q}{m,n} );
                        disp ( hamilt.pol{p,q}{m,n} );
                        grid ( hamilt.pol{p,q}{m,n} );
                    end
                end
            end
        end
    end
else
    prt.disp('***************************************************************')
    prt.disp('No polarizabilities available                     ')
    prt.disp('***************************************************************')
    prt.disp('   ')
end

%% Truncating potential and kinetic energy
hamilt = trunc_pot_kin(hamilt.truncate, hamilt, space)
