%-------------------------------------------------------------------------
%
% Solve the time-independent Schroedinger equation to get eigenstates
% and energies in pseudospectral representation using DVR/FBR techniques 
% 
% Part 3/3: Diagonalizing the Hamiltonian
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2023 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function diag (obj)

% Find eigenvalues and eigenvectors of Hamiltonian matrix
tim = cputime;

% Full matrix: All eigenvalues and eigenvectors
if obj.storage == 'f'
    prt.disp('Start diagonalizing full matrix ...')
    [ obj.eig_vecs, obj.eig_vals ] = eig ( obj.matrix );

% Sparse matrix: Subset of eigenvalues and eigenvectors
else
    prt.disp('Start diagonalizing sparse matrix ...')
    switch obj.choice
        case ('sr') % Smallest real
            [ obj.eig_vecs, obj.eig_vals ] = ...
                eigs ( obj.matrix, obj.number, 'smallestreal' );
        case ('sm') % Smallest magnitude|absolute
            [ obj.eig_vecs, obj.eig_vals ] = ...
                eigs ( obj.matrix, obj.number, 'smallestabs' );
        case ('sigma')
            [ obj.eig_vecs, obj.eig_vals ] = ...
                eigs ( obj.matrix, obj.number, obj.sigma );
    end
end

% Sort eigenvalues of Hamiltonian matrix
eig_vals = diag(obj.eig_vals);
[~, order] = sort(real(eig_vals));
obj.eig_vals = eig_vals(order);

% Sort eigenvectors of Hamiltonian matrix
sortvecs = zeros(size(obj.eig_vecs));
for k= 1:length(order)
    sortvecs(:, k) = obj.eig_vecs(:, order(k));
end
obj.eig_vecs = sortvecs;

prt.disp (['Finished after [CPU seconds] : ' num2str(cputime-tim)])


% Table of eigenvalues
prt.disp(' ')
prt.disp('Table of eigenvalues')
for ii=obj.start : obj.stop
    energy = obj.eig_vals(ii+1);
    prt.disp([ int2str(ii) ' : ' num2str(math.real(energy))])
end
prt.disp(' ')
