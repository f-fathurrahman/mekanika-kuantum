%--------------------------------------------------------------------------
%
% Propagate objects of class wave (wavefunctions on grids) 
% by one substep of size time.steps.s_delta
%
% Fourth order differencing scheme (symplectic, reversible)
% ---------------------------------------------------------
%
%                               i      ^            i      3 ^3                5
%    psi(t+dt) = psi(t-dt) - 2 ---- dt H psi(t) + ------ dt  H  psi(t) + O ( dt ) 
%                              hbar               3*hbar
%
% Initialization in three steps:
% 
% (1) Propagate one quarter timestep backward by explicit Euler scheme
%
%                              i      ^ 
%    psi(t-dt/4) = psi(t) + ------ dt H psi(t)
%                           4*hbar
%
% (2) Propagate one quarter timestep backward by 2nd order differencing
%
%                            i        ^
%    psi(t-dt/2) = psi(t) + ------ dt H psi(t-dt/4) 
%                           2*hbar
%
% (3) Propagate one half timestep backward by 4th order differencing
%
%                          i      ^                  i      3 ^3
%    psi(t-dt) = psi(t) + ---- dt H psi(t-dt/2) - ------- dt  H  psi(t-dt/2)
%                         hbar                    24*hbar
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2023 Burkhard Schmidt's group
%
% see the README file for license details.

classdef diff_4 < handle
   
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = diff_4
           
        end
        
        % Display propagator, overloading default disp method
        function disp(~)
            
            prt.disp ('Finite differences (Euler), symmetrized in time')
            prt.disp ('***************************************************************')
            prt.disp (' ')
            prt.disp ( 'Fourth order differencing ' )
                        
        end


        % Initialize propagator, i.e.
        % provide "old" wavefunction psi(t-dt)
        % Note that psi.dvr must not be changed!
        function init (~, psi)
            global time hamilt

            % Keep psi(t) in memory
            psi.old = psi.dvr;

            % (1) From psi(t) to psi(t-tau/4) by explicit Euler scheme
            apply_ham(psi, [], false);  % psi.dvr -> psi.new
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m} = psi.old{m} + 1i * time.steps.s_delta/4 * psi.new{m};
            end

            % (2) From psi(t-tau/4) and psi(t) to psi(t-tau/2) by 2nd order differencing
            apply_ham(psi, [], false);  % psi.dvr -> psi.new
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m} = psi.old{m} + 1i * time.steps.s_delta/2 * psi.new{m};
            end

            % (3) From psi(t-tau/2) and psi(t) to psi(t-tau) by 2nd order differencing
            apply_ham(psi, [], false);  % psi.dvr -> psi.new
            for m = 1:hamilt.coupling.n_eqs
                psi.sum{m} = psi.old{m} + 1i * time.steps.s_delta * psi.new{m};
            end

            % Apply Hamiltonian two more times
            for h=1:2
                psi.dvr = psi.new;
                apply_ham(psi, [], false); % psi.dvr -> psi.new
            end

            % Fourth order contribution
            for m = 1:hamilt.coupling.n_eqs
                psi.sum{m} = psi.sum{m} - 1/24 * 1i * time.steps.s_delta^3 * psi.new{m};
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

            % Apply Hamiltonian two more times
            for h=1:2
                psi.dvr = psi.new;
                apply_ham(psi, efield, false); % psi.dvr -> psi.new
            end

            % Fourth order contribution
            for m = 1:hamilt.coupling.n_eqs
                psi.sum{m} = psi.sum{m} + 1/3 * 1i * time.steps.s_delta^3 * psi.new{m};
            end

            % Get ready for next time step
            psi.old = psi_t;
            psi.dvr = psi.sum;

        end
    end
end

