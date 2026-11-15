function B = jp_target_blocks(p, nu_list, blk)
% Build the right-hand side b{y,z} of eqs. (8)-(9): for every antenna
% element (column of B) the target phase is piecewise constant in frequency,
% equal to nu_i (+2*pi*k_i, added by the caller) over the subcarriers
% assigned to user i.
%
% nu_list : (M+1)-independent, one Naz x Nel (or Na x 1) matrix per user,
%           already including the 2*pi*k offset.
% blk     : 1 x (M+1) user index per subcarrier
%
% B : (M+1) x Nelem

Nu = numel(nu_list);
Nelem = numel(nu_list{1});
B = zeros(p.M+1, Nelem);
for i = 1:Nu
    sel = (blk == i);
    if ~any(sel), continue; end
    B(sel,:) = repmat(reshape(nu_list{i}, 1, Nelem), nnz(sel), 1);
end
end
