%--------------------------------------------------------------
%
% Get a colormap , trying to match default Matlab curve colors
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2023-... Burkhard Schmidt's group
%
% see the README file for license details.

function cmap = get_colormap(Choice,Length)

    switch Choice
        case 1
            cmap = winter(Length);
        case 2
            cmap = autumn(Length);
        case 3
            cmap = hot(Length);
        case 4
            cmap = cool(Length);
        case 5
            cmap = summer(Length);
        otherwise
            prt.error ('No colormap for this choice')
    end

end
