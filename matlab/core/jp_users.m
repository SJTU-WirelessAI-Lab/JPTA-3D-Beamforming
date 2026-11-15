function [az, el] = jp_users(Nu)
% Sec. IV: user i (i = 0..Nu-1) sits at
%   theta_az,i = -60 + i*120/(Nu-1)  [deg]
%   theta_el,i =  90 + i*30 /(Nu-1)  [deg]
% Note the paper's formula is undefined for Nu = 1; we then use broadside
% of the sector (az = 0, el = 90) which is the natural single-user case.

if Nu < 1 || Nu ~= round(Nu)
    error('Nu must be a positive integer');
end
i = 0:(Nu-1);
if Nu == 1
    az = 0;  el = 90;
else
    az = -60 + i*120/(Nu-1);
    el =  90 + i*30 /(Nu-1);
end
end
