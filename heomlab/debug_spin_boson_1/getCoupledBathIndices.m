function coupled_bath_indices = getCoupledBathIndices(coupled_mode_indices,mode_info)


n_couplings = numel(coupled_mode_indices);
coupled_bath_indices = zeros([1,n_couplings]);
n_debye = mode_info.n_debye;
M = mode_info.M;
n_debye_baths = n_debye / (M+1);
n_bo_baths = 0; %n_bo / (M+2);
N_pade = 0;
n_custom_modes = 0;
n_modes_bath = [(M+1)*ones([1,n_debye_baths]),(M+2)*ones([1,n_bo_baths]), N_pade+1, n_custom_modes];
n_modes_bath_cumsum = cumsum(n_modes_bath);
for n = 1:n_couplings
  jk_index = coupled_mode_indices(n);
  coupled_bath_indices(n) = sum(jk_index > (n_modes_bath_cumsum) ) + 1;
end

end