%--------------------------------------------------------------------------
%
% External electric field as a sequence of shaped pulses with envelope f
% and constant or (linearly or quadratically) chirped carrier frequencies 
% with user-specified polarizations in the x/y plane and
% possibly with phase shifts.
%
%            N
% E(time) = sum f (time) * cos ( omega  * time + phase )
%           i=1  i                    i               i
%
% Returning a cell vector of length two. Each entry contains a vector of 
% field values for any given vector (timesteps) of time values.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function efield = eval (timesteps)
global time

% Preallocate for all polarization directions
efield = cell(time.efield.n_polar,1);
for p=1:time.efield.n_polar
    efield{p} = zeros(size(timesteps));
end

% Loop over pulses
for e = 1:length(time.pulse)
    if ~isempty(time.pulse{e})

        % Get envelopes and oscillations via methods of the pulse objects
        shape = envelope  (time.pulse{e},timesteps);
        vibra = oscillate (time.pulse{e},timesteps);

        % Find those timesteps where pulse is "on"
        where = find( abs(timesteps-time.pulse{e}.delay) <= time.pulse{e}.length/2 );

        % Dressed state (Floquet) picture: Return half envelope
        if time.efield.dressed
            for p=1:time.efield.n_polar
                efield{p}(where) = efield{p}(where) + ...
                    0.5 * shape(where) * time.pulse{e}.ampli(p);
            end

        % Bare state picture: Return oscillating field
        else
            for p=1:time.efield.n_polar
                efield{p}(where) = efield{p}(where) + ...
                    shape(where) .* vibra(where) * time.pulse{e}.ampli(p);
            end
        end

    end

end
end
