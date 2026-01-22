function j_bath = getBathIndexFromModeIndex(jk_index,mode_info)
n_debye = mode_info.n_debye;
M = mode_info.M;
debye_index_max = n_debye;
n_debye_baths = n_debye / (M+1);
if (jk_index <= debye_index_max)
  j_bath = ceil(jk_index/(M+1));
else
  j_bath = n_debye_baths + ceil((jk_index - n_debye)/(M+2));
end
end