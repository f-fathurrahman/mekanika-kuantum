% Solve the time-independent Schroedinger equation to get eigenstates
% and energies in position representation by using DVR/FBR techniques 
%
% Part 1/3: Set up the Hamiltonian matrix in DVR space

function setup(obj)

prt.disp('ENTER ham.eigen.setup')

global hamilt space

prt.disp (' ')
prt.disp('Start setting up matrix ...')
tim = cputime;

% Total number of grid points
N = space.n_tot;

% Initialize matrix: Full or sparse
if obj.storage == 'f'
    obj.matrix = zeros(obj.mat_size);
else
    obj.matrix = sparse(obj.mat_size, obj.mat_size);
end

% Create kinetic energy matrix from both the grid's internal kinetic energy
% operators and external kinetic energy operators. If we solve the TISE for
% coupled equations, the resulting matrix should have a blockwise-diagonal form;
% we fill only the left upper square and copy the result to the other blocks.

% Grid's internal kinetic energy operators
% Optionally already with cutoff and/or sparsification 
for k = 1:space.n_dim
    prt.disp(['... kinetic energy for dof : ' num2str(k)])
    obj.matrix(1:N, 1:N) = ...
        obj.matrix(1:N,1:N) + kinetic2dvr(space.dof{k}, ...
        obj.cutoff, obj.storage);
end

% If available, get external kinetic energy operators
% Cutoff and/or sparsification not yet(!) available 
if isfield(hamilt, 'kin')
    for ii = 1:length(hamilt.kin)
        prt.disp(['... custom kinetic energy operator : ' num2str(ii)]);
        obj.matrix(1:N, 1:N) = ...
            obj.matrix(1:N,1:N) + kinetic2dvr(hamilt.kin{ii});
    end
end

% Copy kinetic energy from upper left squares to other diagonal blocks
for m = 2:hamilt.coupling.n_eqs
    obj.matrix((m-1)*N+1:m*N, (m-1)*N+1:m*N) = ...
        obj.matrix(1:N, 1:N);
end

% Add the potential energy to the diagonal elements of the kinetic energy
prt.disp('... potential energy')
for m = 1:hamilt.coupling.n_eqs
    if ~hamilt.pot{m,m}.empty
        
        % Get potential energy and sparsify if desired
        pot = hamilt.pot{m,m}.dvr(:);
        if obj.storage == 's'
            pot = sparse(pot);
        end
        
        % Add potential energy to the diagonal blocks
        obj.matrix((m-1)*N+1:m*N, (m-1)*N+1:m*N) = ...
            obj.matrix((m-1)*N+1:m*N, (m-1)*N+1:m*N) + diag(pot);
    end
end

% If applicable, get diabatic coupling elements
if hamilt.coupling.n_eqs>1
    prt.disp('... diabatic potential coupling')
    for m = 1:hamilt.coupling.n_eqs
        for n = m+1:hamilt.coupling.n_eqs
            if ~hamilt.pot{m,n}.empty

                % Get potential coupling and sparsify if desired
                pot = hamilt.pot{m,n}.dvr(:);
                if obj.storage == 's'
                    pot = sparse(pot);
                end

                % Diabatic couplings should be symmetric
                obj.matrix((m-1)*N+1:m*N, (n-1)*N+1:n*N) = diag(pot);
                obj.matrix((n-1)*N+1:n*N, (m-1)*N+1:m*N) = diag(pot);

            end
        end
    end
end

% Cut-off and density of matrix (in percent)
obj.matrix(abs(obj.matrix) < obj.cutoff) = 0;
density = 100*nnz(obj.matrix) / obj.mat_size^2;

prt.disp(['Density of matrix : ' num2str(density) ' %'])
prt.disp(['Finished after [CPU seconds] : ' num2str(cputime-tim)])
prt.disp(' ')

prt.disp('EXIT ham.eigen.setup')