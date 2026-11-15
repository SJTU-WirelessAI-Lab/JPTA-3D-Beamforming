function [blk, edges] = jp_block(p, alpha)
% Bandwidth partition of C16 Sec. III-A: user i gets the subcarriers in
% [f_{alpha_{i-1}(M+1)}, f_{alpha_i (M+1)}), i.e. a contiguous block of
% floor/round(alpha_i*(M+1)) subcarriers.
%
% blk   : 1 x (M+1) vector, blk(k) = user index owning subcarrier k
% edges : 1 x (Nu+1) block boundaries (edges(1)=0, edges(end)=M+1)

alpha = alpha(:).';
if any(alpha < 0) || abs(sum(alpha)-1) > 1e-9
    error('alpha must be non-negative and sum to 1 (got %s)', mat2str(alpha));
end
Nu = numel(alpha);
cum = [0, cumsum(alpha)]*(p.M+1);
edges = round(cum);
edges(1) = 0;
edges(end) = p.M+1;
% guard against zero-length blocks
if any(diff(edges) < 1)
    warning('jp_block:emptyBlock', 'a user got 0 subcarriers');
end

blk = zeros(1, p.M+1);
for i = 1:Nu
    if edges(i+1) > edges(i)
        blk(edges(i)+1 : edges(i+1)) = i;
    end
end
end
