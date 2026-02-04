%% SUBROUTINE 1: EXAPND.M
function B = my_expand(A,S)
% Compact version of: http://www.mathworks.co.uk/matlabcentral/fileexchange/24536-expand

% Get the size (and number of dimensions) of input.
SA = size(A);
if length(SA)~=length(S);
    error('Length of size vector must equal ndims(A).  See help.')
elseif any(S~=floor(S));
    error('The size vector must contain integers only.  See help.')
end
for ii = length(SA):-1:1
    H = zeros(SA(ii)*S(ii),1);    % One index vector into A for each dim.
    H(1:S(ii):SA(ii)*S(ii)) = 1;  % Put ones in correct places.
    T{ii} = cumsum(H);            % Cumsumming creates the correct order.
end
B = A(T{:});                      % Feed the indices into A.
