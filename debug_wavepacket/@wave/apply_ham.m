%---------------------------------------------------------------------
%
% Application of Hamiltonian operators to wavefunction objects
%
% Vector efield contains components of external electric field
% or half its envelope (if dressed states are used).
% Argument  norm==1 to use normalized Hamiltonian (see Chebychev)
% 
% Original wavefunction is assumed in property obj.dvr (unchanged!)
% Resulting wavefunction is made available through property obj.new
%
% Input argument "efield" is the electric field 
% (scalar or row vector) which may also be empty
%
% Input argument "normalize" toggles the use of a
% normalized Hamiltonian (Chebychev propagators only)
%
%---------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011      Ulf Lorenz
%
% see the README file for license details.

function apply_ham (obj, efield, normalize)
global  hamilt space

%% Kinetic energy operators
kinetic_grid = cell(hamilt.coupling.n_eqs, 1);
for m = 1:hamilt.coupling.n_eqs
    kinetic_grid{m} = zeros(size(obj.dvr{1}));
end

% Inherent kinetic energy operator for each degree of freedom.
for k = 1:space.n_dim
    kinetic(space.dof{k}, obj, false);
    for m = 1:hamilt.coupling.n_eqs
        kinetic_grid{m} = kinetic_grid{m} + obj.new{m};
    end
end

% Additional kinetic operators
if isfield(hamilt, 'kin')
    for k = 1:length(hamilt.kin)
        kinetic(hamilt.kin{k}, obj, false);
        for m = 1:hamilt.coupling.n_eqs
            kinetic_grid{m} = kinetic_grid{m} + obj.new{m};
        end
    end
end

% Save result of kinetic operators
obj.new = kinetic_grid;
    

%% Potential energy
for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        if ~hamilt.pot{m,n}.empty
            obj.new{m} = obj.new{m} +      hamilt.pot{m,n}.dvr  .* obj.dvr{n};
            if m~=n % Hamiltonian should be Hermitian (symmetric)
                obj.new{n} = obj.new{n} + conj(hamilt.pot{m,n}.dvr) .* obj.dvr{m};
            end
        end
    end
end

%% Dipole moments along all polarization directions
if isfield(hamilt, 'dip') && ~isempty(efield)
    for p = 1:length(hamilt.dip)
        if abs(efield(p))>0
            if ~isempty (hamilt.dip{p})
                for m = 1:hamilt.coupling.n_eqs
                    for n = m:hamilt.coupling.n_eqs
                        if ~hamilt.dip{p}{m,n}.empty
                            obj.new{m} = obj.new{m} -      efield(p) * hamilt.dip{p}{m,n}.dvr  .* obj.dvr{n};
                            if m~=n % Hamiltonian should be Hermitian (symmetric)
                                obj.new{n} = obj.new{n} - conj(efield(p) * hamilt.dip{p}{m,n}.dvr) .* obj.dvr{m};
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Polarizabilities along all polarization directions
if isfield(hamilt, 'pol') && ~isempty(efield)
    for p = 1:size(hamilt.pol,1)
        if abs(efield(p))>0
            for q = p:size(hamilt.pol,2)
                if abs(efield(q))>0
                    if ~isempty (hamilt.pol{p,q})
                        for m = 1:hamilt.coupling.n_eqs
                            for n = m:hamilt.coupling.n_eqs
                                if ~hamilt.pol{p,q}{m,n}.empty
                                    obj.new{m} = obj.new{m} -      efield(p)*efield(q)/2 * hamilt.pol{p,q}{m,n}.dvr  .* obj.dvr{n};
                                    if m~=n % Hamiltonian should be Hermitian (symmetric)
                                        obj.new{n} = obj.new{n} - conj(efield(p)*efield(q)/2 * hamilt.pol{p,q}{m,n}.dvr) .* obj.dvr{m};
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% If desired, use normalized Hamiltonian (for Chebychev propagators only)
if normalize
    for m = 1:hamilt.coupling.n_eqs
        obj.new{m} = obj.new{m} - obj.dvr{m} * ( hamilt.range.delta/2 + hamilt.range.e_min );
        obj.new{m} = obj.new{m} * 2 / hamilt.range.delta;
    end
end
