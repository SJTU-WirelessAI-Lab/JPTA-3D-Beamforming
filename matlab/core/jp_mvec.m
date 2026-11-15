function m = jp_mvec(p)
% Subcarrier index vector, symmetric around DC: -M/2 .. M/2  (M+1 values).
% f_m = fc + m*df.
m = (-p.M/2):(p.M/2);
end
