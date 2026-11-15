function G = jp_gain_freqs(p, phi, tau, az_deg, el_deg, fvec)
% Beamforming gain (linear) of C16 eq. (3) over a vector of frequencies.
%
%   G(f_m) = | sum_{y,z} exp( j*( phi_yz + 2*pi*f_m*tau_yz - Omega_yz(f_m) ) ) |^2 / (Naz*Nel)
%
% with Omega_yz(f) = pi*(f/fc)*( y*sin(az)*sin(el) + z*cos(el) ).
%
% The exact ratio f/fc is used here (no narrowband approximation), because
% this is the physically evaluated gain.  The *design* target b of eq. (9)
% uses nu(y,z) without the ratio, as printed in the paper.
%
% Returns G : 1 x numel(fvec)  (linear, in [0, Naz*Nel])

nu = jp_nu(p, az_deg, el_deg);          % Naz x Nel
nu = nu(:).';                            % 1 x Nelem
phi = phi(:).';  tau = tau(:).';

fvec = fvec(:).';
ratio = fvec / p.fc;                     % 1 x Nk

H  = phi.' + 2*pi*(tau.' * fvec);        % Nelem x Nk
% NOTE: jp_nu already carries the factor pi (nu = pi*(y sin a sin e + z cos e)),
% so Omega = (f/fc)*nu.  Do NOT multiply by pi again here.
Om = nu.' * ratio;                       % Nelem x Nk

AF = sum( exp(1j*(H - Om)), 1 );         % 1 x Nk
G  = abs(AF).^2 / p.Nelem;
end
