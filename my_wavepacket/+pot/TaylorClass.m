%--------------------------------------------------------------------------
% Represent the potential energy function as a Taylor series
% Includes free particle (V=0) as a special case
%--------------------------------------------------------------------------

classdef TaylorClass < pot.generic & handle
    
    properties (Access = public)
        
        hshift      % Horizontal shift (row vector)
        vshift      % Vertical shift (scalar)
        coeffs      % Coefficients, i.e. derivatives
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = TaylorClass()
            obj.empty   = false;
            obj.hshift  = [];
            obj.vshift  = 0;
            obj.coeffs  = [];
        end
        
        % Initialize potential: Set/check parameters
        function init(obj, space)
            if isempty(obj.hshift)
                obj.hshift = zeros(1, space.n_dim);
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            disp@pot.generic (obj)
            prt.disp ('Taylor series (diag. in N dimensions)')
            prt.disp ('***************************************************************')
            prt.disp (' ')
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)  
            V = math.taylor (...
                r, ...
                obj.hshift, ...
                obj.vshift, ...
                obj.coeffs, 1 );
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            F = math.taylor_d (...
                r, ...
                obj.hshift, ...
                obj.coeffs );
        end
        
    end
end
