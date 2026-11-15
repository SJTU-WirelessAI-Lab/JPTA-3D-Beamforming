function [x1, x2] = jp_ls_fit(m, B)
% Exact least-squares solution of A*x = B with A = [1, m] (eq. (7)).
% x = (A'A)^(-1) A'B.  Because the subcarrier index vector m is symmetric
% (M+1 odd, sum(m)=0), A'A = diag(M+1, sum(m.^2)) and the solution is
%   x1 = mean(B)                 (the paper's eq. (11))
%   x2 = sum(m .* B) / sum(m.^2) (slope; the paper's eq. (12),(13))
% x2 is the total phase increment per subcarrier, i.e. 2*pi*df*tau.
%
% B : (M+1) x P  (one column per antenna element / per array axis)
% x1: 1 x P, x2: 1 x P

m = m(:);
assert(numel(m) == size(B,1), 'dimension mismatch');
assert(abs(sum(m)) < 1e-9*max(1,numel(m)), 'm must be symmetric around 0');

x1 = mean(B, 1);
x2 = (m.' * B) / sum(m.^2);
end
