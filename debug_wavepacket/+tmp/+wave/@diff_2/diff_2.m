%--------------------------------------------------------------------------
%
% Propagate objects of class wave (wavefunctions on grids) 
% by one substep of size time.steps.s_delta
%
% Second order differencing scheme (symplectic, reversible)
% ---------------------------------------------------------
%
%                               i      ^                3
%    psi(t+dt) = psi(t-dt) - 2 ---- dt H psi(t) + O ( dt ) 
%                              hbar              
%
%    A. Askar, A. S. Cakmak, J. Chem. Phys. 68, 2794 (1978)
%    https://doi.org/10.1063/1.436072
%
% Initialization in two steps, see page 69 of 
% https://doi.org/10.1016/0021-9991(91)90137-A
% 
% (1) Propagate one half timestep backward by explicit Euler scheme
%
%                              i      ^ 
%    psi(t-dt/2) = psi(t) + ------ dt H psi(t)
%                           2*hbar
%
% (2) Propagate one half timestep backward by 2nd order differencing
%
%                          i      ^  
%    psi(t-dt) = psi(t) + ---- dt H psi(t-dt/2)  
%                         hbar              
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2023 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

classdef diff_2 < handle
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = diff_2
           
        end
        
        % Display propagator, overloading default disp method
        function disp(~)
            
            prt.disp ('Finite differences (Euler), symmetrized in time')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( 'Second order differencing ' )
                        
        end


        % Initialize propagator, i.e.
        % provide "old" wavefunction psi(t-dt)
        % Note that psi.dvr must not be changed!
        function init (~, psi)
            global time hamilt

            % Keep psi(t) in memory
            psi.old = psi.dvr;

            % From psi(t) to psi(t-tau/2) by explicit Euler scheme
            apply_ham(psi, [], false);  % psi.dvr -> psi.new
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m} = psi.old{m} + 1i * time.steps.s_delta/2 * psi.new{m};
            end

            % From psi(t-tau/2) and psi(t) to psi(t-tau) by 2nd order differencing
            apply_ham(psi, [], false);  % psi.dvr -> psi.new
            for m = 1:hamilt.coupling.n_eqs
                psi.sum{m} = psi.old{m} + 1i * time.steps.s_delta * psi.new{m};
            end

            % Get ready for next time step
            psi.dvr = psi.old;
            psi.old = psi.sum;
        end


        % Perform propagation, for given psi(t) and psi(t-tau)
        function propa (~, psi,k)
            global time hamilt

            % Instantaneous electric field
            if isfield(time, 'pulse')
                efield = zeros(1,time.efield.n_polar); % row vector
                for p=1:time.efield.n_polar
                    efield(p) = time.efield.grid{p}(k+time.steps.offset);
                end
            else
                efield = [];
            end

            % Keep psi(t) in memory
            psi_t = psi.dvr;

            % Apply Hamiltonian to wavefunction(s)
            apply_ham(psi, efield, false); % psi.dvr -> psi.new

            % Second order contribution
            for m = 1:hamilt.coupling.n_eqs
                psi.sum{m} = psi.old{m} - 2 * 1i * time.steps.s_delta * psi.new{m};
            end

            % Get ready for next time step
            psi.old = psi_t;
            psi.dvr = psi.sum;

        end
    end
end

