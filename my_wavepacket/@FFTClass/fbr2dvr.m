%--------------------------------------------------------------------------
%
% This routine reconstructs a wavefunction from a plane wave expansion
% using FFTW as included in Matlab) to transform from FBR to DVR
% representation. See grid_dvr2fbr for some more extensive documentation
% of what we do.
%
%--------------------------------------------------------------------------

function dvr = fbr2dvr(obj, fbr)

% Note that we have introduced some prefactors in the expansion, which we have
% to get rid of again.

dvr = fbr ./ obj.kin_factor;
dvr = ifft  ( ifftshift ( dvr, obj.dof ), [], obj.dof );
