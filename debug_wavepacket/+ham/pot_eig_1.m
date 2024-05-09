%--------------------------------------------------------------------------
%
% Get instantaneous (effective) potential energy curve
% for a single channel Schrödinger equation
%
% Input argument "efield" is the electric field 
% (scalar or row vector) which may also be empty
%
% Input argument "save_trafo" toggles saving of 
% the transformation matrices obtained here
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function pot_eig_1 ( efield, save_trafo )
global hamilt space

% Bare potential
if ~hamilt.pot{1,1}.empty 
    dressed = hamilt.pot{1,1}.dvr;
else % free particle
    dressed = zeros (size ( space.dvr{1} ) );
end
    
% Dipole moments along all polarization directions
if isfield(hamilt, 'dip') && ~isempty(efield)
    for p = 1:length(hamilt.dip)
        if abs(efield(p))>0
            if ~isempty (hamilt.dip{p})
                if ~hamilt.dip{p}{1,1}.empty
                    dressed = dressed - efield(p) * hamilt.dip{p}{1,1}.dvr;
                end
            end
        end
    end
end

% Polarizabilities along all polarization directions
if isfield(hamilt, 'pol') && ~isempty(efield)
    for p = 1:size(hamilt.pol,1)
        if abs(efield(p))>0
            for q = p:size(hamilt.pol,2)
                if abs(efield(q))>0
                    if ~isempty (hamilt.pol{p,q})
                        if ~hamilt.pol{p,q}{1,1}.empty
                            dressed = dressed - efield(p)*efield(q)/2 * hamilt.pol{p,q}{1,1}.dvr;
                        end
                    end
                end
            end
        end
    end
end

% Eigenvalues (effective potential)
hamilt.eig_val    = cell (1);
hamilt.eig_val{1} = dressed;

% Transformation matrix: trivial
if save_trafo 
    hamilt.eig_vec    = cell (1);
    hamilt.eig_vec{1,1} = ones(size(dressed));
end

end