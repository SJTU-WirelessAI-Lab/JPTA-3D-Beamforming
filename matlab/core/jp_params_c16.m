function p = jp_params_c16()
% Parameters of C16: Yildiz et al., "3D Beamforming Through Joint Phase-Time
% Arrays", IEEE VTC2024-Fall, pp. 1-7, doi:10.1109/VTC2024-Fall63153.2024.10757955.
% Table I (default values).
%
% All frequencies in Hz, delays in seconds, angles in degrees.

p.tag  = 'C16';
p.fc   = 28e9;          % center frequency
p.BW   = 95e6;          % system bandwidth
p.df   = 120e3;         % subcarrier spacing
p.M    = 792;           % number of subcarriers is M+1 = 793, m = -M/2 .. M/2
p.Naz  = 16;            % antennas in y (azimuth) direction
p.Nel  = 24;            % antennas in z (elevation) direction
p.tau_p   = 2.5e-9;     % delay quantization step
p.tau_max = 200e-9;     % delay range
p.beta_ps = 6;          % phase shifter precision (bits)
p.c    = 299792458;

% Phase-unwrapping convention for the analytical target of eq. (9).
% false : keep the piecewise target continuous across subcarrier blocks by
%         adding per-element 2*pi*K (this is what reproduces Fig.2a);
% true  : use the raw steering value nu, i.e. no integer offsets.
p.nowrap = false;

% derived
p.Nelem = p.Naz * p.Nel;
p.Gmax_dB = 10*log10(p.Nelem);   % 25.84 dB, upper bound used by Algorithm 1
end
