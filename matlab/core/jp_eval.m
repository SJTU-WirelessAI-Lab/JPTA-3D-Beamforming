function R = jp_eval(p, phi, tau, az, el, blk)
% Per-user average beamforming gain of C16 eqs. (4)-(5).
%   Gmean_i = mean_{m in F_i} G(theta_i, f_m)
%   gdb_i   = 10*log10(Gmean_i)
%   Gl      = sum_i gdb_i          (eq. 5, the paper's objective)
%   Gl_avg  = Gl / Nu              (what Fig. 2(a) plots)

m = jp_mvec(p);
f = p.fc + m*p.df;
Nu = numel(az);

R.Gmean  = zeros(1,Nu);
R.gdb    = zeros(1,Nu);
for i = 1:Nu
    sel = (blk == i);
    if ~any(sel), R.Gmean(i) = 1e-12; R.gdb(i) = -120; continue; end
    g = jp_gain_freqs(p, phi, tau, az(i), el(i), f(sel));
    R.Gmean(i) = mean(g);
    R.gdb(i)   = 10*log10(R.Gmean(i));
end
R.Gl     = sum(R.gdb);
R.Gl_avg = R.Gl/Nu;
R.max_gain_over_users = max(R.gdb);
R.user_spread = max(R.gdb) - min(R.gdb);
end
