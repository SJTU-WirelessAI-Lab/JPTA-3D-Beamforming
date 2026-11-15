function run_a1()
% A1 of the repro plan = C16 Fig. 2(a): average user gain vs bandwidth
% allocation ratio for Nu = 2 users, comparing joint / separated designs
% solved with least squares or the infinity norm.
%
% Paper's claims:
%   * joint always better than separated by 0.6 .. 4 dB
%   * mean gap 1.29 dB (LS) and 1.18 dB (infinity norm)
%   * for alpha <= 0.2 the infinity norm beats LS

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));
out = fullfile(here,'out'); if ~exist(out,'dir'), mkdir(out); end

p   = jp_params_c16();
Nu  = 2;
[az, el] = jp_users_c16(Nu);
avec = 0.10:0.05:0.50;
modes = {'joint-inf','joint-ls','sep-inf','sep-ls'};
lbl   = {'Joint Inf Norm','Joint LS [10]','Separated Inf Norm','Separated LS [10]'};

fprintf('===== A1 : C16 Fig.2(a), Nu = %d =====\n', Nu);
fprintf('users at az = [%s] deg, el = [%s] deg\n', num2str(az), num2str(el));
fprintf('alpha grid: 0.10:0.05:0.50 (user 2 gets 1-alpha)\n\n');

G = nan(numel(modes), numel(avec));
T = nan(numel(modes), numel(avec));      % run time per point
span = nan(numel(modes), numel(avec));

for mi = 1:numel(modes)
    for k = 1:numel(avec)
        a  = [avec(k), 1-avec(k)];
        bl = jp_block(p, a);
        t0 = tic;
        [phi, tau, info] = jp_design(p, az, el, a, modes{mi});
        R  = jp_eval(p, phi, tau, az, el, bl);
        T(mi,k) = toc(t0);
        G(mi,k) = R.Gl_avg;              % what Fig. 2(a) plots
        span(mi,k) = info.tau_span_raw*1e9;
    end
end

% ---------- table ----------
fprintf('  alpha |  JointInf  JointLS  SepInf   SepLS  | J-S(LS) J-S(inf)\n');
fprintf('  ------+--------------------------------------+-----------------\n');
for k = 1:numel(avec)
    fprintf('   %.2f |  %7.3f  %7.3f  %7.3f  %7.3f | %+7.2f %+8.2f\n', avec(k), ...
        G(1,k), G(2,k), G(3,k), G(4,k), G(2,k)-G(4,k), G(1,k)-G(3,k));
end

dLS  = G(2,:) - G(4,:);
dINF = G(1,:) - G(3,:);
fprintf('\n  joint - separated (LS ) : min %+.2f  max %+.2f  mean %+.2f dB   (paper: 0.6..4, mean 1.29)\n', min(dLS), max(dLS), mean(dLS));
fprintf('  joint - separated (inf) : min %+.2f  max %+.2f  mean %+.2f dB   (paper: 0.6..4, mean 1.18)\n', min(dINF), max(dINF), mean(dINF));

% fairness check for the alpha <= 0.2 region: inf-norm should beat LS
lo = avec <= 0.2 + 1e-12;
fprintf('\n  alpha <= 0.2 : inf-joint - ls-joint = %s dB\n', mat2str(round(G(1,lo)-G(2,lo),2)));
fprintf('  alpha <= 0.2 : inf-sep   - ls-sep   = %s dB\n', mat2str(round(G(3,lo)-G(4,lo),2)));
fprintf('  mean runtime per design point: joint-ls %.2f s, joint-inf %.2f s, sep-ls %.2f s, sep-inf %.2f s\n', ...
        mean(T(2,:)), mean(T(1,:)), mean(T(4,:)), mean(T(3,:)));
fprintf('  delay span (ns): %s\n', mat2str(round(span(2,:),1)));

% ---------- figure (Fig. 2a) ----------
fg = figure('Visible','off','Position',[100 100 640 480]);
hold on; box on; grid on;
st = {'-o','-s','--^','--d'};
cl = [0.85 0.33 0.10; 0 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56];
for mi = [1 2 3 4]
    plot(avec, G(mi,:), st{mi}, 'Color', cl(mi,:), 'LineWidth',1.6, ...
         'MarkerSize',6, 'MarkerFaceColor', cl(mi,:));
end
xlabel('Bandwidth allocation ratio \alpha (user 1)');
ylabel('Average gain of the users (dB)');
title(sprintf('C16 Fig.2(a) reproduction: N_u = %d, %d\\times%d JPTA', Nu, p.Naz, p.Nel));
legend(lbl, 'Location','southeast');
ylim([18 26]); xlim([0.08 0.52]);
print(fg, fullfile(out,'a1_gain_vs_alpha.png'), '-dpng', '-r180');
close(fg);

% ---------- csv ----------
fid = fopen(fullfile(out,'a1_gain_vs_alpha.csv'),'w');
fprintf(fid,'alpha,joint_inf,joint_ls,sep_inf,sep_ls,tspan_jointls_ns\n');
for k = 1:numel(avec)
    fprintf(fid,'%.2f,%.4f,%.4f,%.4f,%.4f,%.2f\n', avec(k), G(1,k),G(2,k),G(3,k),G(4,k), span(2,k));
end
fclose(fid);
save(fullfile(out,'a1_gain_vs_alpha.mat'),'avec','G','T','span','az','el','p');
fprintf('\nwrote %s\n', fullfile(out,'a1_gain_vs_alpha.png'));
end
