%--------------------------------------------------------------------------
%
% Visualize multidimensional wavepacket or trajectory propagations
% by means of two-dimensional reduced densities which are plotted 
% in DVR representation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-20xy Burkhard Schmidt's group
%
% see the README file for license details.

classdef reduced_2d < vis.densities & handle
    
    properties (Access = public)

        cnt_nlev    % Number of contours: density
        rho_max;    % Maximal values of (coupled) densities 
        
    end
                
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = reduced_2d

            obj = obj@vis.densities;
            
            obj.cnt_nlev  = 20;          % Number of contours for reduced densities
            obj.rho_max;                 % Maximal values of (coupled) densities 
            
        end
        
        %-------------------------------------------------
        % Show reduced densities
        %-------------------------------------------------
        function show_plot (obj, state, step)
            
            global expect space hamilt
            
            if ~isa(state,'wave') && ~isa(state,'traj')
                prt.error('Reduced density (2D) plots available only for wavefunctions or trajectory bundles')
            end
            
            % Even number of dimensions only
            if ~mod(space.n_dim,2)==0
                prt.error ('This visualization for even number of dimensions only')
            end
            
            %% Full Q/M only
            if isa(state,'wave')
                                
                % Calculate reduced densities
                psi.redu = cell (hamilt.coupling.n_eqs, space.n_dim/2);
                
                % Loop over (coupled) wavefunctions
                for m = 1:hamilt.coupling.n_eqs
                    
                    % Only if population exceeds a certain threshold
                    if expect.pop.cha{m}(step)>expect.min_pop
                        
                        for k = 1:space.n_dim/2
                            psi.redu{m,k} = zeros(space.dof{2*k-1}.n_pts,space.dof{2*k}.n_pts);
                        end
                        
                        switch space.n_dim
                            case 2
                                wave_redu_2d(obj,state,m);
                            case 4
                                wave_redu_4d(obj,state,m);
                            case 6
                                wave_redu_6d(obj,state,m);
                           otherwise
                                prt.error ('Reduced Q/M densities only up to 6 dimensions')
                        end
                        
                    end
                    
                end
                
            end
            
            %% Visualize reduced densities
            p = 1;
            m = ceil(sqrt(space.n_dim/2));
            for k = 1:space.n_dim/2
                
                % Arranging the subplots
                if obj.wide % Wide format: 16:9
                    subplot (m,2*m,p)
                    if mod(p,m)==0
                        p = p+1+m;
                    else
                        p = p+1;
                    end
                else % Square format: 9:9
                    subplot (m,m,p)
                    p = p+1;
                end
                if ~obj.hold; hold off; end
                plot( 1234567890, 1234567890 );
                hold on;
                
                % Making the individual plots
                switch(lower(obj.represent))
                    case 'dvr'
                        show_2d_dvr ( obj, state, step, k );
                    otherwise
                        prt.error ('Only the DVR representation is available')
                end
                
            end
            
        end
    end
    
    methods ( Access = private )
        
        show_2d_dvr ( obj, state, step, k );
        
        wave_redu_2d ( obj, state, m)
        wave_redu_4d ( obj, state, m)
        wave_redu_6d ( obj, state, m)
        
    end
end

