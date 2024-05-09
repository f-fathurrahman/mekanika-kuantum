%-------------------------------------------------------------------------
%
% Transform (matrix-valued) potential and (vector-valued) wavefunction
% from diabatic to adiabatic representation or back. It is important 
% that the dia2adi transformation has to be called BEFORE the adi2dia
% because the latter one uses (diabatic) data stored by the former one.
%
% Note that within the context of fully quantum-mechanical WavePacket 
% simulations the adiabatic picture is used only for the purpose
% of visualizing thee wavefunctions and calculating corresponding 
% expectation values but NOT for the propagation itself! This is
% because of the problems associated with the (near) singularities 
% of the kinetic (derivative) couplings at (avoided) crossings and 
% intersections.
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function adiabatic ( psi, step, direction )
global hamilt space time

% Use this method only *during* propagation
if step < 0
    return
end

% Use this method only if adiabatic representation is desired
if strcmpi(hamilt.coupling.represent,'dia')
    return
end

%% Diabatic to adiabatic ("forward transformation")
if strcmpi(direction,'dia2adi')
    
    % Instantaneous electric field
    if isfield(time, 'pulse')
        efield = zeros(1,time.efield.n_polar); % row vector
        for p=1:time.efield.n_polar
            efield(1,p) = time.efield.grid{p}(time.steps.offset);
        end
    else
        efield = [];
    end
    
    % Solve eigenproblem: For time-dependent Hamiltonians
    % for every time step; else only for the first step
    if step==1 || isfield(time, 'pulse')
        
        switch hamilt.coupling.n_eqs
            case 1 % Light-dressed potential
                ham.pot_eig_1 (efield, true)
                
            case 2 % Analytic diagonalization for two coupled channels
                ham.pot_eig_2 (efield, true)
                
            otherwise % Numerical diagonalization for more than two coupled channels
                ham.pot_eig_N (efield, true)
        end
        
    end
    
    % Save diabatic potentials (diagonal only). No necessity to save 
    % diabatic potential couplings because we don't overwrite them
    for m = 1:hamilt.coupling.n_eqs % diag and off-diag
        if ~hamilt.pot{m,m}.empty
            hamilt.pot{m,m}.dia = hamilt.pot{m,m}.dvr;
        end
        for n = m:hamilt.coupling.n_eqs % diag and off-diag
            hamilt.pot{m,n}.dia_empty = hamilt.pot{m,n}.empty;
        end
    end
    
    % Obtain adiabatic potentials: diagonal only, don't overwrite off-diagonals
    for m = 1:hamilt.coupling.n_eqs
        % Note that the class of pot may not be correct any more
        hamilt.pot{m,m}.dvr = hamilt.eig_val{m}; 
        hamilt.pot{m,m}.empty = false; % Adiabatic potentials are never empty
        for n = m+1:hamilt.coupling.n_eqs
            % Adiabatic potential are set to empty (but in reality they are not)
            hamilt.pot{m,n}.empty = true;
        end
    end
    
    % Save diabatic wave functions 
    for m = 1:hamilt.coupling.n_eqs
        psi.dia{m} = psi.dvr{m};
    end
    
    % Transform to adiabatic representation using Hermitian conjugate of eigenvector matrix 
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m} = zeros (size(space.dvr{1}));
        for n = 1:hamilt.coupling.n_eqs
            psi.dvr{m} = psi.dvr{m} + conj(hamilt.eig_vec{n,m}).*psi.dia{n};
        end
    end

    %% Adiabatic to diabatic ("back transformation")
elseif strcmpi(direction,'adi2dia')

    % Retrieve diabatic wavefunctions and potentials (diagonal only)
    for m = 1:hamilt.coupling.n_eqs % diag only
        if ~hamilt.pot{m,m}.dia_empty
            hamilt.pot{m,m}.dvr = hamilt.pot{m,m}.dia; 
        end
        for n = m:hamilt.coupling.n_eqs % diag and off-diag
            hamilt.pot{m,n}.empty = hamilt.pot{m,n}.dia_empty;
        end
    end

    % Retrieve diabatic wavefunctions
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr{m}    = psi.dia{m};
    end

end
