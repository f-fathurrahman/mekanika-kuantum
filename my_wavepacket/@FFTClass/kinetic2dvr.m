% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function dvrkin = kinetic2dvr(obj,cutoff,storage)

prt.disp('')
prt.disp('<<<<< ENTER dof.fft.kinetic2dvr')

global space

% If kinetic operator is disabled
if obj.nokin
	dvrkin = zeros(space.n_tot,space.n_tot);
	return;
end

% Use FGH method to analytically calculate the matrix elements
% See Tannor's book for details.
if obj.periodic
    T_nn = obj.p_max^2 / (6*obj.mass) * (1 + 2/obj.n_pts^2);
else
    T_nn = obj.p_max^2 / (6 * obj.mass);
end
kinetic = eye(obj.n_pts) * T_nn;

for n=1:obj.n_pts
    for m=1:n-1
        if obj.periodic
            T_nm = obj.p_max^2 * (-1)^(n-m) / (obj.mass * obj.n_pts^2 ...
                   * (sin((n-m)*pi/obj.n_pts))^2);
        else
            T_nm = obj.p_max^2*(-1)^(n-m) / ((pi*(n-m))^2 * obj.mass);
        end
        kinetic(n,m) = T_nm;
        kinetic(m,n) = T_nm;
    end
end

% Optionally perform cut-off and convert to sparse 
kinetic( abs(kinetic) < cutoff ) = 0;
if storage == 's'
    kinetic = sparse(kinetic);
    dvrkin = sparse(space.n_tot, space.n_tot);
else
    dvrkin = zeros(space.n_tot, space.n_tot);
end

% Now we basically have to upgrade our matrix T_ij to the generic form
% T_{i1,i2,...j1,j2,...} = eye(i1,i2) * eye(i2,j2) * ... * T_{ik,jk) * ...
% where the asterisk denotes the direct product. Unfortunately, MatLab
% is very poor with multidimensional arrays, and especially seems not to
% have anything like a direct product. Furthermore, if we want to allow
% sparse matrices, Matlab only handles them as two-dimensional objects, which
% prevents an "elegant" reshaping method to be used.
%
% As a basic solution, we setup the 2D matrix, then we go through each
% row, and fill it.

indices = cell(space.n_dim, space.n_dim);
[indices{:}] = ind2sub(size(space.dvr{1}), (1:space.n_tot));

for ii = 1:space.n_tot
    allowed    = ones(space.n_tot, 1);

    % Now disallow all indices, where the other dimension's indices change,
    % which translates to a zero matrix element
    for k = 1:space.n_dim
        if k == obj.dof
            continue;
        end

        allowed(indices{k} ~= indices{k}(ii)) = 0;
    end

    % get the indices of the allowed cells, and fill the row
    allowedindex = find(allowed);
    dvrkin(ii, allowedindex) = kinetic( indices{obj.dof}(ii), indices{obj.dof}(allowedindex) );
end

disp('<<<<< EXIT dof.fft.kinetic2dvr')