classdef FFTClass < handle
    
    % ffr: all properties are made public

    properties (Access = public)
        
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass that enters the kinetic energy
        periodic    % use periodic boundary conditions or not
        nokin       % enable/disable the kinetic energy operator
        
        n_pts       % number of grid points (equal in position and momentum grid)
        x_min       % lower bound of the position grid
        x_max       % upper bound of the position grid
        
        weight      % weights in DVR ("position space")
        x_grid      % grid points in DVR ("position space")
        p_grid      % grid points in FBR ("momentum space")

        x_dlt       % (constant!) grid spacing in position space
        p_dlt       % (constant!) grid spacing in momentum space
        
        p_min       % lower bound of the momentum grid
        p_max       % upper bound of the momentum grid
        
        kin         % DVR grid representation of the kinetic energy operator
        kin_shift   % an fftshifted version of kin
        kin_expo    % same as grid, but exponentiated, for split operator
        kin_factor  % normalize and removes a spurious complex phase
        
        kin_max
    end
    
    methods (Access = public) % in separate files within *this* directory
        
        dvr = fbr2dvr(obj, fbr)

        fbr = dvr2fbr(obj, dvr)
        
        init_grid(obj)
        
        %init_kin(obj, fraction, output)
        init_kin(obj, space, time_var, hamilt, fraction, output)
        
        kinetic(obj, psi, new)
        
        kinetic_exp(obj, psi)
        
        %dvrkin = kinetic2dvr(obj,cutoff,storage)
        dvrkin = kinetic2dvr(obj, space, cutoff, storage)
        
        dvr = matrix2dvr(obj, fbr)
        
        retval = momentum(obj, psi)

        % ffr: These methods are disabled
        % obj = subsasgn(obj, index, val)
        % retval = subsref(obj, s)
        
        % Constructor method: Initialize default property values
        function obj = FFTClass()
            obj.dof = 1;
            obj.mass = 1;
            obj.periodic = true;
            obj.nokin = false;
        end
        
        function disp(obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            prt.disp('***************************************************************')
            prt.disp([ 'DVR for the degree of freedom: ' obj.label])
            prt.disp('Discretization scheme: Equally spaced grids (FFT)')
            prt.disp('***************************************************************')
            prt.disp(' ')
            prt.disp([ 'Number of grid points  : ' num2str(obj.n_pts) ])
            prt.disp(' ')
            prt.disp('Position space')
            prt.disp([ 'Minimum of grid values : ' num2str(obj.x_min) ])
            prt.disp([ 'Maximum of grid values : ' num2str(obj.x_max) ])
            prt.disp([ 'Spacing of grid values : ' num2str(obj.x_dlt) ])
            prt.disp(' ')
            prt.disp('Momentum (wavenumber) space (Fourier transform)')
            prt.disp([ 'Maximum of grid values : ' num2str(obj.p_max) ])
            prt.disp([ 'Spacing of grid values : ' num2str(obj.p_dlt) ])
            prt.disp(' ')
            prt.disp('Momentum (wavenumber) space (Wigner transform)')
            prt.disp([ 'Maximum of grid values : ' num2str(obj.p_max/2) ])
            prt.disp([ 'Spacing of grid values : ' num2str(obj.p_dlt/2) ])
            prt.disp(' ')
            
        end
        
        % Evaluate kinetic energy function (for trajectories!)
        function eval_kin(obj, state)
            state.kin = state.kin + state.mom{obj.dof}.^2 / (2*obj.mass);
        end
        
    end
    
end
