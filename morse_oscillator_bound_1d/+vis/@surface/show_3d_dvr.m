%--------------------------------------------------------------
%
% Visualize wavepacket in 3 dimensions in position (DVR) space 
% using Matlab's builtin iso-surface plots
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2023 Burkhard Schmidt's group
%
% see the README file for license details.

function show_3d_dvr ( obj, state, step )
global expect hamilt info space

if isa(state,'traj')
    prt.error ('Code missing for surface plots of trajectory data in 3D')
end
    
% Deletes all visible graphics objects from previous time steps
if step>1
    cla 
end

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs

    % Optionally plot potential energy surfaces
    if obj.energy && ~hamilt.pot{m,m}.empty

        % Set iso-energy values (linearly spaced)
        if obj.srf_color(2)
            iso_values = linspace(obj.col_min(2), obj.col_max(2), obj.srf_nlev(2));
        else
            iso_values = linspace(hamilt.pot_min, hamilt.pot_max/2, obj.srf_nlev(2));
        end

        % Optionally use a color map; else use default colors
        if obj.srf_maps(2)
            cmap = obj.get_colormap(m,length(iso_values));
        end

        % Loop over iso-surface levels
        for iso=1:obj.srf_nlev(2)
            pot = patch(isosurface ( space.dvr{1},space.dvr{2},space.dvr{3}, ...
                hamilt.pot{m,m}.dvr, iso_values(iso)));
            if obj.srf_maps(2)
                set(pot, ...
                    'FaceColor', cmap(iso,:), ... 
                    'EdgeColor','none', ...
                    'FaceAlpha',0.25);   % High transparency
            else
                set(pot, ...
                    'FaceColor', obj.colors(m,:), ...
                    'EdgeColor','none', ...
                    'FaceAlpha',0.25);   % High transparency
            end

            % Use equal data unit lengths along each axis
            axis equal 
        end

    end
    
    % Plot densities from wavefunctions
    if expect.pop.cha{m}(step)>expect.min_pop

        % Get position density
        rho = abs ( state.dvr{m} ) .^2;

        % Set iso-density values (linearly spaced)
        if obj.srf_color(1)
            iso_values = linspace(obj.col_min(1), obj.col_max(1), obj.srf_nlev(1));
        else
            iso_values = linspace(obj.scale_dvr*0.1, obj.scale_dvr*0.9, obj.srf_nlev(1));
        end

        % Optionally use a color map; else use default colors
        if obj.srf_maps(1)
            cmap = obj.get_colormap(m,length(iso_values));
        end

        % Loop over iso-surface levels
        for iso=1:obj.srf_nlev(1)
            p = patch(isosurface ( space.dvr{1},space.dvr{2},space.dvr{3}, ...
                rho, iso_values(iso)));
            if obj.srf_maps(1)
                set(p, ...
                    'FaceColor', cmap(iso,:), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.75)       % Low transparency
            else
                set(p, ...
                    'FaceColor', obj.colors(m,:), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.75)       % Low transparency
            end

            % Use equal data unit lengths along each axis
            axis equal
        end
    end
    
    if m==1
        hold on
    end
    
end
hold off

% Specify view point in terms of azimuth and elevation
view (obj.srf_view(1),obj.srf_view(2))

% Lighting of surfaces
if strcmpi(info.system,'Matlab')
    camlight
    lighting gouraud
    material shiny
elseif strcmpi(info.system,'Octave')
    camlight;
end
if obj.srf_look(2)
    if strcmpi(info.system,'Matlab')
        lightangle(obj.srf_light(1),obj.srf_light(2));
    elseif strcmpi(info.system,'Octave')
        camlight(obj.srf_light(1),obj.srf_light(2));
end

% Axes, labels, etc
if ~obj.range
    axis ([ ...
        space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{2}.dvr_min space.dof{2}.dvr_max ...
        space.dof{3}.dvr_min space.dof{3}.dvr_max ]);
else
    axis ([ ...
        obj.x_min obj.x_max ...
        obj.y_min obj.y_max ...
        obj.z_min obj.z_max ]);
end

set(gca, ...
    'LineWidth',  obj.l_thick, ...
    'FontName',   obj.f_name, ...
    'FontSize',   obj.f_large, ...
    'FontWeight', obj.f_heavy)

title  ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
zlabel ( [ 'R_{', space.dof{3}.label, '}' ] )

end
