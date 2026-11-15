function nu = jp_nu(p, az_deg, el_deg)
% Array steering value for the (y,z)-th element, eq. (9):
%   nu(y,z,i) = y*pi*sin(az_i)*sin(el_i) + z*pi*cos(el_i)
% Returned as an Naz x Nel matrix.

y = (0:p.Naz-1).';
z = (0:p.Nel-1);
az = az_deg*pi/180;
el = el_deg*pi/180;

nu = pi*( y * (sin(az)*sin(el)) ) + pi*( z * cos(el) );
end
