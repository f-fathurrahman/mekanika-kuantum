%--------------------------------------------------------------------
% Input of a bilinear control system: external electric field vs. time
%--------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2019-20.. Burkhard Schmidt
%
% see the README file for license details.

function input (obj,~)
global time

if time.efield.dressed
    factor = 2;
else
    factor = 1;
end

% Set range of time (sub!)steps
if time.steps.offset == 1
    steps = 1; % first main time step
else
    steps = 1:time.steps.offset+time.steps.s_number; % all later main time steps
end

% Plot electric field(s)
if isfield(time,'pulse')
    for p=1:time.efield.n_polar % all polarization directions
        h = plot ( time.steps.s_grid(steps), factor * real(time.efield.grid{p}(steps)) );
        set(h, ...
            'LineStyle',   obj.patterns{p}, ...
            'Color',       'black', ...
            'LineWidth',   obj.l_thick, ...
            'DisplayName', ['F_' int2str(p)])
        if p==1
            hold on
        end
    end
    hold off
end

% Legend explaining the line styles
if obj.legends 
    legend('Location','SouthWest')
end

% Axes and labels
if time.efield.dressed
    axis ( [ 0 time.steps.t_total  [-1.1 +1.1]*time.efield.max_ampli ] )
else
%     line ([0 time.steps.t_total],[0 0],'Color','black','LineWidth', obj.l_thin,'DisplayName','F=0')
    axis ( [ 0 time.steps.t_total -abs(time.efield.max_ampli)*1.1 +abs(time.efield.max_ampli)*1.1 ] )
end

% Fonts and line styles
set(gca, ...
	'LineWidth',  obj.l_thick, ...
	'FontName',   obj.f_name, ...
	'FontSize',   obj.f_large, ...
	'FontWeight', obj.f_heavy, ...
	'Box',        'on')

% Axis labels
xlabel ('t' )
ylabel ('input: fields' )

end
        
