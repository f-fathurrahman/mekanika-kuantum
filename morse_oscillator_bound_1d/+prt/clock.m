% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function clock
global info

% CPU time elapsed
tt = cputime-info.start_time;

% Date and time (nicely formatted)
dd = date;
c = clock;
hh = num2str(c(4),'%2.2i');
mm = num2str(c(5),'%2.2i');
ss = num2str(c(6),'%5.2f');

% Output
prt.disp ('***************************************************************');
if tt<0.1
    prt.disp ( 'Starting the clock ...' )
else
    prt.disp (['Elapsed seconds : ' num2str(tt,'%11.5e') ' seconds'] );
end
prt.disp ('***************************************************************');
prt.disp ( ' ' )
prt.disp (['Date            : ', dd] )
prt.disp (['Time            : ', hh, ':', mm, ':', ss] );
prt.disp ( ' ' )
