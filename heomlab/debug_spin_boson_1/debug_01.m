clear variables; close all;

% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters
epsilon = 1.0;
Delta = 2.0;
% bath parameters
beta = 1.0;
% debye bath parameters
lambda_D = 0.5;
omega_D = 2.0;

% dynamics information - parameters for the Short-Iterative Arnoldi
% Integrator
dt = 1e-2;
n_steps = 1000;
krylov_dim = 8;
krylov_tol = 1e-8;
% parmeters for heirarchy truncation using L/M truncation
L_max = 6; 
M_max = 4;


% matrices of system observable operators to be returned, sigma_x, sigma_y
% sigma_z, and 1
O_sys = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)};

% initial state of the system
rho_0_sys = [[1,0];[0,0]];

% two objects are supplied to the HEOM dynamics function:
% "full_system" specifies the full Hamiltonian (system + bath) and the
% temperature.
% "heom_dynamics" specifies the HEOM truncation, integrator for the
% dynamics and the total propagation time and observables to be calculated.

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon,Delta];
                     [Delta,-epsilon]];
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",[[1,0];[0,-1]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)};
full_system.beta = beta;

% a struct that contains information about the HEOM dynamics
heom_dynamics = struct;
% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct;
heom_dynamics.integrator.method = "SIA";
heom_dynamics.integrator.dt = dt;
heom_dynamics.integrator.n_steps = n_steps;
heom_dynamics.integrator.krylov_dim = krylov_dim;
heom_dynamics.integrator.krylov_tol = krylov_tol;

% hierarchy trunction information
heom_dynamics.heom_truncation = struct;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off";
heom_dynamics.heom_truncation.M_max = M_max;
heom_dynamics.heom_truncation.L_max = L_max;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off";
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut;
heom_dynamics.heom_truncation.heom_termination = "markovian";

% what system observables should be returned
heom_dynamics.observables = struct;
heom_dynamics.observables.system = O_sys;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys;

% run the dynamics
%[O_t,t] = runHEOMDynamics(full_system, heom_dynamics);

% ...Call to runHEOMDynamics will be "inlined" here


% convert the bath info into a more use-able form
%heom_bath_info = getBathInformation(full_system);
% ffr: there are other variables that should be returned but are discarded?

% getBathInfomation
fprintf('\n---- ENTER getBathInformation ----\n');

% ffr: why also return lambda_Ds, etc ?
baths = full_system.baths;
n_baths = numel(baths);

lambda_Ds = [];
omega_Ds = [];
lambda_OBOs = [];
Omega_OBOs = [];
gamma_OBOs = [];
lambda_UBOs = [];
Omega_UBOs = [];
gamma_UBOs = [];
lambda_Ds_pade = [];
omega_Ds_pade = [];
pade_approximants = [];
N_pade = [];
Vs = {};
nus_custom = {};
cs_custom = {};
cbars_custom = {};
cs_trunc_custom = {};
nus_trunc_custom = {};
cbars_trunc_custom = {};

for i = 1:n_baths
  fprintf('Pass here 119\n')
  if (baths{i}.spectral_density == "debye")
    lambda_Ds = [lambda_Ds,baths{i}.lambda_D];
    omega_Ds = [omega_Ds,baths{i}.omega_D];
    Vs = [Vs,{baths{i}.V}]; 
  end
end

heom_bath_info = struct();
heom_bath_info.n_baths = n_baths;
heom_bath_info.Vs = Vs;
heom_bath_info.lambda_Ds = lambda_Ds;
heom_bath_info.omega_Ds = omega_Ds;
heom_bath_info.lambda_OBOs = lambda_OBOs;
heom_bath_info.Omega_OBOs = Omega_OBOs;
heom_bath_info.gamma_OBOs = gamma_OBOs;
heom_bath_info.lambda_UBOs = lambda_UBOs;
heom_bath_info.Omega_UBOs = Omega_UBOs;
heom_bath_info.gamma_UBOs = gamma_UBOs;
heom_bath_info.beta = full_system.beta;
heom_bath_info.lambda_Ds_pade = lambda_Ds_pade;
heom_bath_info.omega_Ds_pade = omega_Ds_pade;
heom_bath_info.N_pade = N_pade;
heom_bath_info.pade_approximants = pade_approximants;
heom_bath_info.nus_custom = nus_custom;
heom_bath_info.cs_custom = cs_custom;
heom_bath_info.cbars_custom = cbars_custom;
heom_bath_info.nus_trunc_custom = nus_trunc_custom;
heom_bath_info.cs_trunc_custom = cs_trunc_custom;
heom_bath_info.cbars_trunc_custom = cbars_trunc_custom;

fprintf('\n---- EXIT getBathInformation ----\n');



fprintf('\n---- ENTER constructHEOMGenerator ----\n');

% construct the HEOM dynamics genrator as a sparse matrix
%[L_heom, ado_indices] = constructHEOMGenerator(full_system.H_sys, heom_bath_info, ...
%    heom_dynamics.heom_truncation);

% There are 3 arguments for this case
H_sys = full_system.H_sys;
% heom_bath_info
heom_truncation_info = heom_dynamics.heom_truncation;

% set up the hierarchy using the Djkstra frequency cut-off method using
% matsubara decompositions of the bath correlation functions

% This is the case of heom_truncation_info.truncation_method == "depth cut-off")
% get the specified maximum matsurbara mode
M = heom_truncation_info.M_max;
% get the max depth of the hierarchy
L_max = heom_truncation_info.L_max;


% get arrays the ado frequencies (nus) and coupling coefficents (cs)
% for the debye, OBO and UBO baths
nus = [];
cs = [];
cbars = [];
lambdas = [];
mode_info = struct();
mode_info.M = M;

% This is the case of lambda_Ds of bath_info are given
n_debye = numel(heom_bath_info.lambda_Ds);
[nus_array_debye,cs_array_debye] = generateNusAndCsDebye(heom_bath_info.omega_Ds,...
    heom_bath_info.lambda_Ds,heom_bath_info.beta,M);
n_debye = numel(nus_array_debye);
nus = [nus,reshape(transpose(nus_array_debye),[1,n_debye])];
cs = [cs,reshape(transpose(cs_array_debye),[1,n_debye])];
cbars = [cbars,reshape(transpose(cs_array_debye),[1,n_debye])];
mode_info.n_debye = n_debye;
lambdas = [lambdas, reshape(repmat(heom_bath_info.lambda_Ds,[M+1,1]),[1,n_debye])];


% No need for OBO for in this case
mode_info.n_obo = 0;
% No need for UBO for in this case
mode_info.n_ubo = 0;
% No need for Debye Pade
n_debye_pade_baths = 0;
mode_info.n_debye_pade = 0;
mode_info.N_pade = [];
cs_array_debye_pade = [];
nus_array_debye_pade = [];
% No need for nu_custom
mode_info.n_custom_baths = 0;
mode_info.n_custom_modes = [];



% This is the case of heom_truncation_info.truncation_method == "depth cut-off")
% construct the hierarchy structure with depth (L) based truncation
[ado_indices,ado_gammas,lower_indices,upper_indices, ...
  coupled_mode_indices,truncated_coupled_modes, ...
  ado_indices_term,modes_term,term_indices] = generateHierarchyDepthTrunc(L_max,nus);

% create an array of the coupled mode indices
coupled_bath_indices = getCoupledBathIndices(coupled_mode_indices,mode_info);
heom_structure = struct();
heom_structure.ado_gammas = ado_gammas;
heom_structure.ado_indices = ado_indices;
heom_structure.lower_indices = lower_indices;
heom_structure.upper_indices = upper_indices;
heom_structure.coupled_mode_indices = coupled_mode_indices;
heom_structure.truncated_coupled_modes = truncated_coupled_modes;
heom_structure.coupled_bath_indices = coupled_bath_indices;
heom_structure.nus = nus;
heom_structure.cs = cs;
heom_structure.cbars = cbars;
heom_structure.M = M;
heom_structure.mode_info = mode_info;
heom_structure.ado_indices_term = ado_indices_term;
heom_structure.modes_term = modes_term;
heom_structure.term_indices = term_indices;
heom_structure.cs_array_debye = cs_array_debye;
heom_structure.nus_array_debye = nus_array_debye;
heom_structure.mode_info = mode_info;
heom_structure.ado_indices_term = ado_indices_term;
heom_structure.modes_term = modes_term;
heom_structure.term_indices = term_indices;


% get some dimensions of things for setting up the HEOM generator
d_sys = size(H_sys,1);
d_liou = d_sys^2;
n_ados = size(ado_indices,1);
d_heom = d_liou*n_ados; % total dimension of the HEOM system of ADOs
V = heom_bath_info.Vs;
beta = heom_bath_info.beta;
n_baths = length(V);
n_debye_baths = numel(heom_bath_info.lambda_Ds);
n_OBO_baths = numel(heom_bath_info.lambda_OBOs);
n_UBO_baths = numel(heom_bath_info.lambda_UBOs);
n_debye_pade_baths = numel(heom_bath_info.lambda_Ds_pade);
n_couplings = size(lower_indices,1);

fprintf('N_ado = %d, M = %d\n',[n_ados,M]);

% first make some useful system superoperators
id_sys = speye(d_sys);
id_liou = speye(d_liou);
id_ados = speye(n_ados);

% free system liouvillian
L_sys = -1.0i*(kron(H_sys,id_sys) - kron(id_sys,transpose(H_sys)));

% superoperators that correspond to left (L) and right (R) multiplication
% by the bath operators and the renormalisation operator
V_L = {};
V_R = {};
V_comm = {};
for j = 1:n_baths
  V_L{j} = kron(V{j},id_sys);
  V_R{j} = kron(id_sys,transpose(V{j}));
  V_comm{j} = V_L{j}-V_R{j};
end

% add matsurbara truncation correction
% This is the case of (heom_truncation_info.heom_termination == "markovian"  || heom_truncation_info.heom_termination == "low temp correction")
Xi = sparse([],[],[],d_liou,d_liou);

for j = 1:n_debye_baths
  R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - ...
    heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) - ...
    sum(cs_array_debye(j,2:end)./nus_array_debye(j,2:end));
  Xi = Xi - R_j * V_comm{j}*V_comm{j};
end
Xi = kron(id_ados,Xi);

% add the Pade white-noise term (skipped)

% add the free system evolution term to the HEOM generator
L_heom = kron(id_ados,L_sys) + Xi;

% add the decay terms -sum_jk n_jk nu_jk
L_heom = L_heom - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou);

heom_structure.Xi = Xi;
heom_structure.L_sys = L_sys;
heom_structure.V = V;
heom_structure.V_comm = V_comm;
heom_structure.V_L = V_L;
heom_structure.V_R = V_R;
heom_structure.d_heom = d_heom;
heom_structure.n_ados = n_ados;
heom_structure.d_liou = d_liou;

% construct the HEOM generator
for r = 1:n_couplings
  J = lower_indices(r);
  K = upper_indices(r);
  jk_coup = coupled_mode_indices(r);
  % need to fix getting the bath index!!!!
  %     j_coup = ceil(jk_coup/(M+1)); % get the bath index that is coupling J & K
  j_coup = coupled_bath_indices(r);
  
  J_block = ((J-1)*(d_liou)+1):(J*d_liou);
  K_block = ((K-1)*(d_liou)+1):(K*d_liou);
  n_jk_coup = ado_indices(K,jk_coup);
  c_jk_coup = cs(jk_coup);
  cbar_jk_coup = cbars(jk_coup);
  % add the coupling from the ADO deeper in the heirarchy
  % d/dt rho_J = ... -i [V_j , rho_K]
  L_heom(J_block,K_block) = L_heom(J_block,K_block)...
      -1.0i *sqrt(n_jk_coup*abs(c_jk_coup))* V_comm{j_coup};
  % add the coupling from the ADO shallower in the heirarchy
  if (abs(c_jk_coup)>0)
      L_heom(K_block,J_block) = L_heom(K_block,J_block) ...
          + (-1.0i*c_jk_coup*sqrt(n_jk_coup/abs(c_jk_coup))) * V_L{j_coup} ...
          + (1.0i*conj(cbar_jk_coup)*sqrt(n_jk_coup/abs(c_jk_coup))) * V_R{j_coup};
  end
end

% ORIGINAL INCORRECT TERMINATOR CODE - diagonal only
diag_only_term = false;
% diag_only_term handling is removed

% This is the case of heom_truncation_info.heom_termination == "markovian" || heom_truncation_info.heom_termination == "NZ2")
% add a perturbative correction for the ados at which the hierarchy is terminated
%     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes);
%     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info);
%     n_term = size(terminator_ado_indices,1);
% get the number of terminating ADOs
n_term_ados = size(ado_indices_term,1);
for r = 1:n_term_ados
  % get the number of connections up into the explicit hierarchy
  n_term_indices = numel(term_indices{r});
  for ind_J = 1:n_term_indices
    ind_K_range = 1:n_term_indices;
    for ind_K = ind_K_range
      % get the hierarchy indices of the coupled terms in the
      % truncated hierarchy
      J = term_indices{r}(ind_J);
      K = term_indices{r}(ind_K);
      % get the block indices
      J_block = ((J-1)*(d_liou)+1):(J*d_liou);
      K_block = ((K-1)*(d_liou)+1):(K*d_liou);
      % get the terminating modes
      term_mode_J = modes_term{r}(ind_J);
      term_mode_K = modes_term{r}(ind_K);
      % get terminating bath indices
      j_J = getBathIndexFromModeIndex(term_mode_J,mode_info);
      j_K = getBathIndexFromModeIndex(term_mode_K,mode_info);
      % the the coupling coefficients
      c_J = cs(term_mode_J);
      c_K = cs(term_mode_K);
      cbar_K = cbars(term_mode_K);
      % get the mode excitation levels
      n_J = ado_indices(J,term_mode_J);
      n_K = ado_indices(K,term_mode_K);
      gamma_ado_term = sum(nus.*ado_indices_term(r,:));
      % This is the case of (heom_truncation_info.heom_termination == "markovian")
      L_heom(J_block,K_block) = L_heom(J_block,K_block) - ...
        (sqrt((n_J+1)*abs(c_J))/gamma_ado_term) * V_comm{j_J}*...
          (c_K*sqrt((n_K+1)/abs(c_K)) * V_L{j_K} - ...
           conj(cbar_K)*sqrt((n_K+1)/abs(c_K)) * V_R{j_K} );
    end
  end
end

fprintf('\n---- EXIT constructHEOMGenerator ----')

% Output: 
% L_heom <- L_heom
% ado_indices <- heom_structure
ado_indices = heom_structure; % ???

% get the dimensions of things
d_heom = size(L_heom,1);
d_hilb = size(full_system.H_sys,1);
d_liou = d_hilb * d_hilb;

% construct the rho_0 for the full hierarchy
rho_0_heom = zeros([d_heom,1]);
rho_0_heom(1:d_liou) = convertToLiouvilleVector(heom_dynamics.rho_0_sys);
if isfield(heom_dynamics,"rho_0_heom")
  rho_0_heom = heom_dynamics.rho_0_heom;
end

% set up the observable operators
n_obs = numel(heom_dynamics.observables.system);
O = sparse([],[],[],n_obs,d_heom);
for n = 1:n_obs
  O(n,1:d_liou) = convertToLiouvilleVector(heom_dynamics.observables.system{n})';
end

integrator = heom_dynamics.integrator;
[O_t,t,rho_final_heom] = runDynamicsSIADensityOperator( ...
  rho_0_heom,L_heom,...
  integrator.n_steps, integrator.dt, O, integrator.krylov_dim, ...
  integrator.krylov_tol);

% Outputs: O_t, t, rho_final_heom, L_heom

% plot sigma_alpha(t)
%ylabels = {'\langle\sigma_x(\itt\rm)\rangle','\langle\sigma_y(\itt\rm)\rangle','\langle\sigma_z(\itt\rm)\rangle'};
%figure
%for i = 1:3
%  subplot(3,1,i)
%  plot(t,O_t(i,:))
%  xlabel('\itt\rm')
%  ylabel(ylabels{i})
%end