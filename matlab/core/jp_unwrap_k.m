function K = jp_unwrap_k(p, az, el)
% Wrapping-offset integer k(y,z,i) of eq. (10):
%   k(y,z,1) = 0
%   k(y,z,i) = k(y,z,i-1) + round( (nu_i-1 - nu_i)/(2*pi) )
% Purpose: keep the phase gap between consecutive users below 2*pi.
%
% K is Naz x Nel x Nu.
%
% A per-axis variant (used by the "separated" design) is obtained by calling
% this function with a 1x1 array (Nel = 1 or Naz = 1).

Nu = numel(az);
K = zeros(p.Naz, p.Nel, Nu);
if Nu == 1
    return;
end

nu = zeros(p.Naz, p.Nel, Nu);
for i = 1:Nu
    nu(:,:,i) = jp_nu(p, az(i), el(i));
end
for i = 2:Nu
    K(:,:,i) = K(:,:,i-1) + round( (nu(:,:,i-1) - nu(:,:,i))/(2*pi) );
end
end
