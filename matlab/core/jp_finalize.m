function [phiq, tauq, info] = jp_finalize(p, x1, x2, shape)
% Turn an LS / inf-norm solution (x1 = phase at m=0, x2 = phase slope per
% subcarrier) into physical, quantized PS phases and TTD delays.
%
% The hardware model of C16 eq. (2) is  omega = exp(j(phi + 2*pi*f_m*tau)).
% Writing f_m = fc + m*df gives
%       phi + 2*pi*f_m*tau = (phi + 2*pi*fc*tau) + 2*pi*df*tau*m
% so the LS unknowns map to
%       tau      = x2 / (2*pi*df)
%       phi_PS   = x1 - 2*pi*fc*tau          (PS absorbs the carrier phase)
% Note: NOT phi = x1.  Using phi = x1 would add a spurious term
% (fc/df)*x2 ~ 2.3e5*x2 rad, i.e. a completely different beam.
%
% A common delay offset tau0 added to all elements multiplies the array
% factor by exp(j*2*pi*(f_m-fc)*tau0), i.e. a common phase - the gain is
% unchanged.  We therefore shift tau to be non-negative (hardware needs
% tau >= 0) and record whether the delay range is exceeded.
%
% Outputs are quantized: tau on the tau_p grid in [0, tau_max],
% phi on the 2*pi/2^beta grid.

tau = x2 / (2*pi*p.df);

% shift to non-negative
tau = tau - min(tau(:));

info.tau_span_raw = max(tau(:));                 % required delay range
info.saturated    = max(tau(:)) > p.tau_max + eps;
tau = min(tau, p.tau_max);

% quantize delay first, then set the PS phase so that the total phase at
% m = 0 is preserved
tauq = round(tau / p.tau_p) * p.tau_p;
tauq = min(max(tauq, 0), p.tau_max);

phi = x1 - 2*pi*p.fc*tauq;

q = 2*pi / 2^p.beta_ps;
phiq = mod(round(phi/q)*q, 2*pi);

phiq = reshape(phiq, shape);
tauq = reshape(tauq, shape);
end
