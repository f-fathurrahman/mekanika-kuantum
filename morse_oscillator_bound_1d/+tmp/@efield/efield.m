%------------------------------------------------------------------------------
%
% Temporal discretization: Dealing with electric fields
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2010-2011 Ulf Lorenz
%
% see the README file for license details.

classdef efield < handle
    
    properties (Access=public)
        
        n_polar     % number of polarization directions
        n_pulse     % number of E-field pulses

        dressed     % toggle deressed vs. bare states
        photons     % dressed with this many photons
        complex     % toggle use of complex-valued fields
        
        max_ampli   % Maximum amplitude
        grid        % values of el. field for all time (sub!) steps
        
    end
    
    methods (Access=public)
        
        % Constructor: Set trivial default values
        function obj = efield
            
            obj.complex = false;         % default: not using complex fields
            obj.dressed = false;         % default: not using dressed states
            
        end
        
        % Initialization: Check / set parameters
        function init (obj)
            
            global hamilt time
            
            % Error checking
            if obj.complex && hamilt.coupling.n_eqs > 2
                prt.error('Cannot use complex electric fields with more than 2 coupled states');
            end
            
            % Number of electric field pulses
            obj.n_pulse = length(time.pulse);

            % Number of polarization directions
            obj.n_polar = [];
            for e=obj.n_pulse
                if ~isempty(time.pulse{e})
                    if isempty(obj.n_polar)
                        obj.n_polar = length(time.pulse{e}.ampli);
                    else
                        if obj.n_polar ~= length(time.pulse{e}.ampli)
                            prt.error ('All E-field amplitude vectors should be of the same length')
                        end
                    end
                end
            end
            
            % Maximal amplitude (all E-field pulses, all polarization directions)
            obj.max_ampli = 0;
            for e=1:obj.n_pulse
                if ~isempty(time.pulse{e})
                    for p=1:obj.n_polar
                        if abs(time.pulse{e}.ampli(p)) > abs(obj.max_ampli)
                            obj.max_ampli = time.pulse{e}.ampli(p);
                        end
                    end
                end
            end

            % Set up grid for values of the electric field (all polarizations)
            for p=1:obj.n_polar
                obj.grid{p} = zeros(length(time.steps.s_grid),1);
            end
            
            % Evaluate electric field for all time (sub!) steps
            obj.grid = efi.eval (time.steps.s_grid);
            
        end
        
        % Floquet function: see extra file
        floquet (obj)
        
    end
    
end
