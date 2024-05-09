%--------------------------------------------------------------------------
%
% Generic properties of all potential energy class definitions
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef generic < handle
    
    properties (Access = public)
        row         % row index (in diabatic potential matrix)
        col         % column index (in diabatic potential matrix)
        dvr         % Grid representation (in N dimensions)
        dia         % Backup copy for dia=>adi transformation
        empty       % Is this potential empty?
        dia_empty   % Is this potential empty?
    end
    
    methods (Access = public)
        
        % Constructor: Set empty potential as the default
        function obj = generic
            obj.empty = true;
            obj.dia_empty = true;
        end

        % Initialize: dummy method
        function init(~)
        end

        % Partial copy of an instance of this class, at least a few properties
        % This is used in WavePacket in +tmp/@efield/floquet.m
        function obj2 = copy(obj1)
            obj2 = pot.generic(); % constructor of THIS class
            obj2.row = obj1.row;
            obj2.col = obj1.col;
            obj2.dvr = obj1.dvr;
            obj2.empty = obj1.empty;
        end
        
        % Display potential, overloading default disp method
        function disp (obj)
            global hamilt
            prt.disp('***************************************************************')
            if obj.row == obj.col
                prt.disp(['Potential energy for channel ' hamilt.coupling.labels{obj.row} ':'])
            else
                prt.disp(['Potential coupling for channels ' hamilt.coupling.labels{obj.row} ', ' hamilt.coupling.labels{obj.col} ':'])
            end
            if obj.empty
                prt.disp ('Not available')
                prt.disp ('***************************************************************')
                prt.disp (' ')
            end
        end
        
        % Grid representation of potential energy function
        function grid (obj)
            global space
            
            if ~obj.empty
                obj.dvr = V ( obj, space.dvr );
            end
            
        end
        
    end
    
end

