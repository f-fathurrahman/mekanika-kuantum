% Wave functions represented on DVR / FBR grids
% for use in fully quantum-mechanical dynamics

classdef WaveClass < generic & handle
    
    properties (Access = public)
        
        dvr         % Discrete variable representation (cell vector for coupled states)
        adi         % Backup copy for adi=>dia transformation
        dia         % Backup copy for dia=>adi transformation
        fbr         % Finite basis representation (cell vector for coupled states)
        ini         % Initial wavefunction (DVR, cell vector)
        new         % New wavefunction (DVR, cell vector)
        old         % Old wavefunction (DVR, cell vector)
        sum         % Sum of wavefunctions (DVR, Chebychev only, cell vector)
        mom         % Momentum operator acting on wavefunction (DVR, single channel)
        
        bound       % Bound states; only for relaxation with 'cheby_imag'
        
        redu        % Reduced densities (for plots only)
        wig         % Wigner transform (for plots only, cell vector)
        wig_max     % Maximum of Wigner transform
        
        M_ham       % Matrix (eigen) representation: Hamiltonian 
        M_amo       % Matrix (eigen) representation: Additional mult. op.'s
        M_dip       % Matrix (eigen) representation: Dipole moment
        M_pol       % Matrix (eigen) representation: Polaritability
        M_sbc       % Matrix (eigen) representation: System-bath coupling
        M_mat       % Matrix (eigen) representation: Observables
        M_vec       % Vector (eigen) representation: Observables
        M_lab       % Labels of observables
        M_obs       % Types of observables
        
    end
    
    methods (Access = public)
        
        % Constructor: Setting default values
        function obj = WaveClass()
            % Inherit from constructor of generic superclass
            obj = obj@generic;
            % ffr: why need generic class?
            % ffr: probably it will be used in generic save/load utilities

            % Setting the name of this class
            obj.string0 = 'wave';
            % ffr: need to be renamed to WaveClass ?
            % ffr: Probably not, this is used to differentiate between different kind of representation
            
            % Cell array with bound states; only for 'cheby_imag' 
            obj.bound = {};
            %
        end
 
        % Logfile/console output
        function disp(obj)
            prt.disp('This is a WaveClass object')
        end
        
        % How much memory needed to save a "wave" object
        function out = memory_size(obj, hamilt, space)
            out = numel(space.dvr{1}) * 16 * hamilt.coupling.n_eqs;
        end
        
        % More methods: see separate files
        % ffr: If the method signatures are changed then it also
        % should be changed here

        init_obj(obj)                % Initial conditions

        propagate(obj, step)          % Propagation
        
        eigen(obj, step)          % Eigenfunctions of Hamiltonian
        
        observe(obj, step)          % Expectation values / uncertainties
        
        [hamilt, space, time_var] = init_ham(obj, hamilt, space, time_var)

        apply_ham(obj, efield, norm)  % Application of Hamiltonian
        
        adiabatic(obj, step, direction) % Adiabatic<=>diabatic transformation
        
        diabatic(obj)                % Adiabatic=>diabatic (initial only)
        
        shift_mom(obj, mom_0)         % Shift in momentum space
        
        save(obj, step)          % Saving dynamics to file(s)
        
        load(obj, step)          % Loading dynamics from file(s)
        
        save_0(obj)          % Saving general info to file
        
        load_0(obj, choice)        % Loading general from file
        
    end
    
    methods (Static)
        
        wig = wigner(dvr)
        
        psi_norm = normalize(psi_in)
        retval = braket(bra, ket)
        retval = sandwich(bra, operator, ket)
        
    end % methods Static
    
end

