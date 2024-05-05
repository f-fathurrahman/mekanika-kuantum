function retval = momentum(obj, psi, space)
% global space

if isempty(obj.kin)
    init_kin(obj, 1);
end

% Applies the momentum operator to the input argument psi.
% Since this is usually done only for expectation values, i.e.
% not very often, we do not care about optimisation too much.

% Transform the wave function to FBR
retval = dvr2fbr(obj, psi);

% Apply the momentum. We assume it is stored in space.fbr*
retval = retval .* space.fbr{obj.dof};

% Transform back
retval = fbr2dvr(obj, retval);
