%--------------------------------------------------------------------------
% Potential energy: Morse oscillator
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2009 Burkhard Schmidt
%
% see the README file for license details.

classdef MorseClass < pot.generic & handle
    
    properties (Access = public)
        
        d_e         % dissociation energy
        r_e         % equilibrium position
        alf         % range parameter
        t_e         % energetic shift
                
    end
    
    methods (Access = public)
        
       % Constructor: Set default values
        function obj = MorseClass()
            disp('ENTER morse constructor')
            obj.empty = false;
            obj.r_e = 0;
            obj.t_e = 0;
            disp('EXIT morse constructor')
        end
        
        
        % Initialize potential: Set/check parameters
        function init(obj, space)
            if length (obj.d_e) ~= space.n_dim
                prt.error ('Wrong length of Morse parameters')
            end
            if length (obj.r_e) ~= space.n_dim
                prt.error ('Wrong length of Morse parameters')
            end
            if length (obj.alf) ~= space.n_dim
                prt.error ('Wrong length of Morse prameters')
            end
            if any (obj.d_e<0) || any (obj.r_e<0) || any (obj.alf<0) 
                prt.error ('Morse parameters should not be negative')
            end
        end

        % This is to avoid overloading disp function
        function disp_info(obj, space)
            prt.disp('Morse oscillator object')
            prt.disp('***************************************************************')
            for k = 1:space.n_dim
                prt.disp(' ')
                if space.n_dim > 1
                    prt.disp( [ 'Dimension : ' num2str(k) ])
                end
                prt.disp( [ 'Dissociation energy    : ' num2str(obj.d_e(k)) ] )
                prt.disp( [ 'Equilibrium position   : ' num2str(obj.r_e(k)) ] )
                prt.disp( [ 'Range parameter (alfa) : ' num2str(obj.alf(k)) ] )
            end
            prt.disp( [ 'Energetic shift        : ' num2str(obj.t_e) ] )
            prt.disp(' ')
            % ffr: Scalar in MATLAB can be indexed (?)
            % ffr: pot1.d_e(1) is the same as pot1.d_e
        end

        % Display potential, overloading default disp method
        function disp(obj)
            prt.disp('Morse oscillator object')
        end
        
        % Evaluate potential energy function
        function V = V(obj, r)
            % ffr: r is assumed to be a cell array ?????
            %
            % Summing up contributions from each component of position vector
            V = zeros(size(r{1}));
            
            for k = 1:length(r)
                V = V + obj.d_e(k) ...
                    * ( 1 - exp ( -obj.alf(k) * ( r{k}-obj.r_e(k) ) ) ).^2 ...
                    + obj.t_e;
            end
            
        end
        
        % Evaluate forces as negative gradients of potential
        % yet to be generalized to the multidimensional case
        function F = F(obj,r)
            F{1} = - 2 * obj.d_e * obj.alf ...
                * ( 1 - exp ( -obj.alf * ( r{1}-obj.r_e ) ) ) ...
                .*      exp ( -obj.alf * ( r{1}-obj.r_e ) )  ;
        end
    end
end

