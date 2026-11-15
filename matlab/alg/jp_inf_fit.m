function [x1, x2] = jp_inf_fit(m, B)
% Infinity-norm fit  min_{x1,x2} || [1 m]*[x1;x2] - b ||_inf  for each column
% of B (C16 Sec. III-A-2).  The paper suggests an LP toolbox; the remote
% compute node has no Optimization Toolbox, so this uses the equivalent
% one-dimensional convex formulation, which is exact and needs only core
% MATLAB:
%
%   for a given slope x2 the best intercept is x1 = (L+U)/2 and the achieved
%   error is t = (L-U)/2, where
%        L(x2) = max_i (b_i - m_i*x2),   U(x2) = min_i (b_i - m_i*x2).
%   t(x2) is convex (L is convex, -U is convex), so a bounded scalar
%   minimiser over x2 returns the minimax solution.
%
% Returns x1, x2 : 1 x P

m = m(:);
P = size(B, 2);
x1 = zeros(1,P);  x2 = zeros(1,P);

lo = -1e3;  hi = 1e3;
opt = optimset('TolX', 1e-13, 'Display','off');

for j = 1:P
    b = B(:,j);
    f = @(s) max(b - m*s) - min(b - m*s);
    s = fminbnd(f, lo, hi, opt);
    L = max(b - m*s);  U = min(b - m*s);
    x2(j) = s;
    x1(j) = (L + U)/2;
end
end
