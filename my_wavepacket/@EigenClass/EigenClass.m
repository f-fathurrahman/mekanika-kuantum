%--------------------------------------------------------------------------
%
% Set up matrix representation of the Hamiltonian operator 
% to numerically solve the corresponding eigenproblem, i.e. 
% the time-independent Schrödinger (TISE) equation using the  
% Fourier Grid Hamiltonian (FGH) method.
%
% Originally discovered by R.Meyer 
%   https://doi.org/10.1063/1.1673259
% Later rediscovered by C.C.Marston and G.G.Balint-Kurti
%   https://doi.org/10.1063/1.456888 
% Generalized to the multidimensional case
% see Eq. (13) in work by S. P. Webb and Sh. Hammes-Schiffer
%   https://doi.org/10.1063/1.1289528
%
%--------------------------------------------------------------------------

classdef EigenClass < handle
    
    properties (Access=public)
        
        start       % Index of first eigenstate wanted
        stop        % Index of last eigenstate wanted
        number      % Number of eigenstates wanted
        storage     % Full or sparse storage of Hamiltonian matrix
        cutoff      % Set matrix entries below this threshold to zero
        choice      % How to work with 'eigs'
        sigma       % Compute eigenvalues near this number
        symmetry    % Symmetry (irr. rep.) label
        matrix      % Matrix representation of Hamiltonian
        mat_size    % Size of Hamiltonian matrix
        transform   % Transformation matrix (if symmetry adation)
        
        eig_vals    % Eigenvalues
        eig_vecs    % Eigenvector
        
    end
    
    methods (Access=public)
        
        % Constructor: Set default values
        function obj = EigenClass
            prt.disp('ENTER hamilt.eigen constructor')
            obj.start   = 00;            % start at ground state
            obj.stop    = 10;            % stop at 10th excited state
            obj.cutoff  = eps;           % cutoff at machine epsilon
            obj.storage = 'f';           % using full matrices
            obj.choice  = 'sr';          % eigenvalues with smallest real
            obj.sigma   = 0;             % eigenvalues near this number
            prt.disp('EXIT hamilt.eigen constructor')
        end
        
        function time_var = init(obj, hamilt, space, time_var)
            % global hamilt space time_var
            
            % Artificial mapping of eigenstates onto time step 
            time_var.steps.m_start  = obj.start;
            time_var.steps.m_stop   = obj.stop;
            obj.number          = obj.stop - obj.start + 1;
            time_var.steps.m_number = obj.number;
            time_var.steps.m_grid   = (obj.start : obj.stop)';
            time_var.steps.t_total  = obj.stop;
            time_var.steps.s_delta  = 1e-10; % dummy setup to avoid crashes

            % Matrix size
            obj.mat_size = space.n_tot * hamilt.coupling.n_eqs;
        end
            
        % Display objects of this class
        function disp(obj)
            
            prt.disp('***************************************************************')
            prt.disp('Solve TISE by direct diagonalization of the Hamiltonian  ')
            prt.disp('***************************************************************')
            prt.disp('  ')
            prt.disp(['Size of matrix                      : ' int2str(obj.mat_size) '^2'])
            if obj.cutoff == 0
                prt.disp('No cutoff value');
            else
                prt.disp(['Cut-off absolute values below       : ' num2str(obj.cutoff)])
            end
            if obj.storage == 'f'
                prt.disp('Working with full matrices (eig)');
            elseif obj.storage == 's'
                 switch obj.choice
                    case ('sr')
                        prt.disp('Working with sparse matrices (eigs) : smallest real')
                    case ('sm')
                        prt.disp('Working with sparse matrices (eigs) : smallest magnitude')
                    case ('sigma')
                        prt.disp(['Working with sparse matrices (eigs) : sigma = ' num2str(obj.sigma)])
                    otherwise
                        prt.error('Wrong choice for eigs options')
                end
            else
                prt.error('Wrong storage scheme. hamilt.eigen.storage has to be "f" (full) or "s" (sparse)');
            end
            prt.disp('   ')
            prt.disp(['Lowest eigenstate investigated  : ' num2str(obj.start)])
            prt.disp(['Highest eigenstate investigated : ' num2str(obj.stop)])
             
        end
        
        % Other functions, defined in different files within the current directory
        %setup(obj)
        setup(obj, hamilt, space)

        symm(obj)
        
        diag(obj)
        
    end
    
end

