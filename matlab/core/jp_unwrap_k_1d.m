function K = jp_unwrap_k_1d(nu_cell)
% Same wrapping-offset recursion as jp_unwrap_k, but for a single array axis
% (used by the "separated" design , Sec. III-B).
% nu_cell{i} is a Na x 1 vector.

Nu = numel(nu_cell);
Na = numel(nu_cell{1});
K = zeros(Na, Nu);
for i = 2:Nu
    K(:,i) = K(:,i-1) + round( (nu_cell{i-1} - nu_cell{i})/(2*pi) );
end
end
