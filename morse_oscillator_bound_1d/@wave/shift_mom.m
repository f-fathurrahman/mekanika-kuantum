%--------------------------------------------------------------------------
%
% Shift (coupled) wavefunction(s) in momentum space
% by multiplication with exp(i*k_0*r) where k_0 is a
% momentum vector to be specified by the user.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2023-.... Burkhard Schmidt
%
% see the README file for license details.

function shift_mom (obj,mom_0)
global space hamilt

if length(mom_0) ~= length(space.dvr)
    prt.error ('Wrong length of momentum vector')
end

for m = 1:hamilt.coupling.n_eqs
    for d = 1:space.n_dim
        obj.dvr{m} = obj.dvr{m} .* exp ( 1i * mom_0(d) * space.dvr{d} );
    end
end

end