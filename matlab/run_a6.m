function run_a6()
% A6 of the repro plan = Fig. 5 (Sec. IV-B): greedy algorithm (Alg. 1)
% and Adam gradient descent (Sec. III-D), both initialised from the LS
% solution, for the joint and the separated parameterisation.
%
% Paper's claims:
%   * Fig. 5(a), Nu = 4, equal BW: joint greedy +2.26 dB, joint GD +2.11 dB
%     (stated as "average gain increase" over the analytical LS solution)
%   * Fig. 5(b), Nu = 2, various alpha: greedy and GD comparable, same
%     joint/separated pattern as Fig. 2
%   * GD runs ~250x faster than greedy in the joint case
%
% Parts:
%   0  incremental vs literal Algorithm 1 (numerical equivalence)
%   1  Fig. 5(a)
%   2  Fig. 5(b)
%   3  runtime / the 250x claim
%   4  why the "+2.26 dB" depends on the analytical baseline

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));
out = fullfile(here,'out');  if ~exist(out,'dir'), mkdir(out); end

p = jp_params();
ZETA_G = 1e-4;               % greedy: Alg.1 line 18, per pass
oG = struct('zeta',ZETA_G,'maxIter',40);
oD = struct('mode','joint','lr',0.1,'zeta',1e-5,'win',20,'maxIter',600);

fprintf('===== A6 : Fig.5, greedy (Alg.1) + Adam GD =====\n');
fprintf('array %dx%d = %d elements, fc = %.0f GHz, BW = %.0f MHz, M+1 = %d\n', ...
        p.Naz, p.Nel, p.Nelem, p.fc/1e9, p.BW/1e6, p.M+1);
fprintf('greedy: tau %.1f ns step / %d pts, phase %d bits / %d pts, zeta %g\n', ...
        p.tau_p*1e9, numel(0:p.tau_p:p.tau_max), p.beta_ps, 2^p.beta_ps, ZETA_G);
fprintf('GD    : Adam lr %.2f (rad), zeta %g over a %d-iteration window\n\n', ...
        oD.lr, oD.zeta, oD.win);

%% ---------------- Part 0 : fast vs literal Algorithm 1 ----------------
ps = p;  ps.Naz = 6;  ps.Nel = 8;  ps.Nelem = ps.Naz*ps.Nel;
NuS = 4;  [azS, elS] = jp_users(NuS);  aS = ones(1,NuS)/NuS;
fprintf('--- Part 0: incremental vs literal Algorithm 1 (array %dx%d, 3 passes) ---\n', ps.Naz, ps.Nel);
chk = struct();
for md = {'joint','sep'}
    m1 = md{1};
    o = struct('mode',m1,'zeta',1e-9,'maxIter',3);
    [~,~,i1] = jp_greedy(ps, azS, elS, aS, o);
    [~,~,i2] = jp_greedy_naive(ps, azS, elS, aS, o);
    chk.(m1) = struct('fast',i1.Gl,'naive',i2.Gl,'diff',abs(i1.Gl-i2.Gl), ...
                      'tfast',i1.time,'tnaive',i2.time,'nIter',i1.nIter);
    fprintf('  %-5s : fast Gl %.10f (%.1f s) | literal Gl %.10f (%.1f s) | diff %.2e\n', ...
            m1, i1.Gl, i1.time, i2.Gl, i2.time, chk.(m1).diff);
end

%% ---------------- Part 1 : Fig. 5(a) ----------------
Nu = 4;  [az, el] = jp_users(Nu);  alpha = ones(1,Nu)/Nu;
blk = jp_block(p, alpha);
fprintf('\n--- Fig.5(a): Nu = %d, equal BW (alpha_i = %.2f) ---\n', Nu, alpha(1));

names = {'LS joint (init)','greedy joint','GD joint','LS sep (init)','greedy sep','GD sep'};
G5a = zeros(numel(names), Nu);  T5a = zeros(numel(names),1);  ITS = zeros(numel(names),1);

t1 = tic;  [ph, ta] = jp_design(p, az, el, alpha, 'joint-ls');
R = jp_eval(p, ph, ta, az, el, blk);  G5a(1,:) = R.gdb;  T5a(1) = toc(t1);
t1 = tic;  [ph, ta] = jp_design(p, az, el, alpha, 'sep-ls');
Rs = jp_eval(p, ph, ta, az, el, blk);  G5a(4,:) = Rs.gdb;  T5a(4) = toc(t1);

o = oG; o.mode = 'joint';
[ph, ta, ig] = jp_greedy(p, az, el, alpha, o);
R = jp_eval(p, ph, ta, az, el, blk);  G5a(2,:) = R.gdb;  T5a(2) = ig.time;  ITS(2) = ig.nIter;

[ph, ta, id] = jp_gd(p, az, el, alpha, oD);
R = jp_eval(p, ph, ta, az, el, blk);  G5a(3,:) = R.gdb;  T5a(3) = id.time;  ITS(3) = id.nIter;

o = oG; o.mode = 'sep';
[ph, ta, igs] = jp_greedy(p, az, el, alpha, o);
R = jp_eval(p, ph, ta, az, el, blk);  G5a(5,:) = R.gdb;  T5a(5) = igs.time;  ITS(5) = igs.nIter;

o = oD; o.mode = 'sep';
[ph, ta, ids] = jp_gd(p, az, el, alpha, o);
R = jp_eval(p, ph, ta, az, el, blk);  G5a(6,:) = R.gdb;  T5a(6) = ids.time;  ITS(6) = ids.nIter;

% values read off the published Fig.5(a) (+-0.2 dB), for direct comparison
pap5a = struct('jointGD',[23.50 21.85 21.70 23.45], 'sepGD',[22.20 21.10 20.15 21.30], ...
               'jointGreedy',[23.35 22.25 21.90 23.05], 'sepGreedy',[21.70 21.10 21.10 22.35]);

fprintf('\n  %-14s | %s |  mean  |  iters |  time\n', 'design', ...
        sprintf('%8s%8s%8s%8s','u1','u2','u3','u4'));
fprintf('  ---------------+--------------------------------+--------+--------+--------\n');
for k = 1:numel(names)
    fprintf('  %-14s | %7.3f %7.3f %7.3f %7.3f  | %6.3f | %5d  | %6.2f s\n', ...
            names{k}, G5a(k,:), mean(G5a(k,:)), ITS(k), T5a(k));
end
fprintf('  ---------------+--------------------------------+--------+--------+--------\n');
fprintf('  paper greedy j | %7.2f %7.2f %7.2f %7.2f  | %6.2f |\n', pap5a.jointGreedy, mean(pap5a.jointGreedy));
fprintf('  paper GD     j | %7.2f %7.2f %7.2f %7.2f  | %6.2f |\n', pap5a.jointGD, mean(pap5a.jointGD));
fprintf('  paper greedy s | %7.2f %7.2f %7.2f %7.2f  | %6.2f |\n', pap5a.sepGreedy, mean(pap5a.sepGreedy));
fprintf('  paper GD     s | %7.2f %7.2f %7.2f %7.2f  | %6.2f |\n', pap5a.sepGD, mean(pap5a.sepGD));

fprintf('\n  max |mine - paper|: greedy joint %.2f dB, GD joint %.2f dB, greedy sep %.2f, GD sep %.2f\n', ...
        max(abs(G5a(2,:)-pap5a.jointGreedy)), max(abs(G5a(3,:)-pap5a.jointGD)), ...
        max(abs(G5a(5,:)-pap5a.sepGreedy)), max(abs(G5a(6,:)-pap5a.sepGD)));

fprintf('\n  gain increase over the LS start (average over users):\n');
fprintf('    joint : greedy %+0.2f dB (paper +2.26) | GD %+0.2f dB (paper +2.11)\n', ...
        mean(G5a(2,:))-mean(G5a(1,:)), mean(G5a(3,:))-mean(G5a(1,:)));
fprintf('    sep   : greedy %+0.2f dB | GD %+0.2f dB\n', ...
        mean(G5a(5,:))-mean(G5a(4,:)), mean(G5a(6,:))-mean(G5a(4,:)));
fprintf('  alternative baselines for the joint delta:\n');
fprintf('    vs sep-LS (%.2f dB): greedy %+0.2f, GD %+0.2f\n', mean(G5a(4,:)), ...
        mean(G5a(2,:))-mean(G5a(4,:)), mean(G5a(3,:))-mean(G5a(4,:)));
fprintf('  GD/joint vs greedy/joint : %+0.2f dB   (paper: GD slightly above)\n', ...
        mean(G5a(3,:))-mean(G5a(2,:)));
fprintf('  map-back loss (Gl_drop)  : joint %.1e, sep %.1e\n', id.Gl_drop, ids.Gl_drop);
fprintf('  delay span: greedy j %.1f ns, GD j %.1f ns, greedy s %.1f ns, GD s %.1f ns\n', ...
        ig.tau_span_raw*1e9, id.tau_span_raw*1e9, igs.tau_span_raw*1e9, ids.tau_span_raw*1e9);

% ---- Fig. 5(a) ----
fg = figure('Visible','off','Position',[100 100 700 470]);
hold on; box on; grid on;
cl = [0.55 0.55 0.55; 0.85 0.33 0.10; 0 0.45 0.74; 0.72 0.72 0.72; 0.47 0.67 0.19; 0.49 0.18 0.56];
st = {'--o','-o','-v','--^','-^','-d'};
% the paper's Fig.5(a) shows the four iterative designs only; the LS starts
% are kept in the tables/CSV
kshow = [2 3 5 6];
for k = kshow
    plot(1:Nu, G5a(k,:), st{k}, 'Color', cl(k,:), 'LineWidth',1.7, ...
         'MarkerSize',7, 'MarkerFaceColor', cl(k,:));
end
xlabel('User index');  ylabel('Average gain of the user (dB)');
title(sprintf('Fig.5(a) reproduction: N_u = %d, equal BW', Nu));
legend(names(kshow), 'Location','southeast');
xlim([0.7 Nu+0.3]);  ylim([18 25.6]);
print(fg, fullfile(out,'a6_fig5a_gain_per_user.png'), '-dpng', '-r180');
close(fg);

%% ---------------- Part 2 : Fig. 5(b) ----------------
Nu2 = 2;  [az2, el2] = jp_users(Nu2);
avec = 0.10:0.05:0.50;
lbl2 = {'LS joint','greedy joint','GD joint','LS sep','greedy sep','GD sep'};
G5b = nan(numel(lbl2), numel(avec));
T5b = nan(numel(lbl2), numel(avec));
ITb = nan(numel(lbl2), numel(avec));

fprintf('\n--- Fig.5(b): Nu = %d, alpha sweep (paper: joint GD ~24.6..24.9, joint greedy ~24.4..24.9) ---\n', Nu2);
fprintf('  alpha | LSjoint GrJoint GDJoint | LSsep  GrSep  GDSep | tGr   tGD  | GD-Gr\n');
for k = 1:numel(avec)
    a = [avec(k), 1-avec(k)];
    bl = jp_block(p, a);
    t1 = tic; [ph, ta] = jp_design(p, az2, el2, a, 'joint-ls');
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(1,k)=R.Gl_avg; T5b(1,k)=toc(t1);
    t1 = tic; [ph, ta] = jp_design(p, az2, el2, a, 'sep-ls');
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(4,k)=R.Gl_avg; T5b(4,k)=toc(t1);

    o = oG; o.mode = 'joint'; o.maxIter = 25;
    [ph, ta, ii] = jp_greedy(p, az2, el2, a, o);
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(2,k)=R.Gl_avg; T5b(2,k)=ii.time; ITb(2,k)=ii.nIter;

    o = oD; o.mode = 'joint';
    [ph, ta, ii] = jp_gd(p, az2, el2, a, o);
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(3,k)=R.Gl_avg; T5b(3,k)=ii.time; ITb(3,k)=ii.nIter;

    o = oG; o.mode = 'sep'; o.maxIter = 25;
    [ph, ta, ii] = jp_greedy(p, az2, el2, a, o);
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(5,k)=R.Gl_avg; T5b(5,k)=ii.time; ITb(5,k)=ii.nIter;

    o = oD; o.mode = 'sep';
    [ph, ta, ii] = jp_gd(p, az2, el2, a, o);
    R = jp_eval(p, ph, ta, az2, el2, bl);  G5b(6,k)=R.Gl_avg; T5b(6,k)=ii.time; ITb(6,k)=ii.nIter;

    fprintf('   %.2f | %6.3f  %7.3f  %7.3f | %6.3f %6.3f  %6.3f | %5.1f %5.1f | %+5.2f\n', ...
            avec(k), G5b(1,k),G5b(2,k),G5b(3,k), G5b(4,k),G5b(5,k),G5b(6,k), ...
            T5b(2,k), T5b(3,k), G5b(3,k)-G5b(2,k));
end
fprintf('\n  mean |GD-greedy| joint = %.3f dB, max %.3f  (paper: "comparable")\n', ...
        mean(abs(G5b(3,:)-G5b(2,:))), max(abs(G5b(3,:)-G5b(2,:))));
fprintf('  mean |GD-greedy| sep   = %.3f dB, max %.3f\n', ...
        mean(abs(G5b(6,:)-G5b(5,:))), max(abs(G5b(6,:)-G5b(5,:))));
fprintf('  greedy - LS (joint): %s dB\n', mat2str(round(G5b(2,:)-G5b(1,:),2)));

% ---- Fig. 5(b) ----
fg = figure('Visible','off','Position',[100 100 700 470]);
hold on; box on; grid on;
for k = 1:numel(lbl2)
    plot(avec, G5b(k,:), st{k}, 'Color', cl(k,:), 'LineWidth',1.6, ...
         'MarkerSize',6, 'MarkerFaceColor', cl(k,:));
end
xlabel('Bandwidth allocation ratio \alpha (user 1)');
ylabel('Average gain of the users (dB)');
title(sprintf('Fig.5(b) reproduction: N_u = %d, various BW allocations', Nu2));
legend(lbl2, 'Location','southeast');
ylim([18 25.6]);
print(fg, fullfile(out,'a6_fig5b_gain_vs_alpha.png'), '-dpng', '-r180');
close(fg);

%% ---------------- Part 3 : runtime ----------------
fprintf('\n--- runtime / the "250x" claim (Nu = 4, Fig.5a design point) ---\n');
phT = ones(p.Naz,p.Nel);  taT = zeros(p.Naz,p.Nel);
tE = 0;
for r = 1:5
    t1 = tic; jp_eval(p, phT, taT, az, el, blk); tE = tE + toc(t1);
end
tE = tE/5;
nCand = numel(0:p.tau_p:p.tau_max) + 2^p.beta_ps;
naiveIter = tE * nCand * p.Nelem;
fprintf('  one eq.(5) evaluation at full size : %.1f ms\n', tE*1e3);
fprintf('  literal Alg.1  : %d gain evals per pass -> %.0f s per pass\n', nCand*p.Nelem, naiveIter);
fprintf('  incremental    : %.2f s per pass (%d passes) -> %.1f s total\n', ...
        ig.time/max(ig.nIter,1), ig.nIter, ig.time);
fprintf('  Adam GD        : %.3f s per iteration (%d iterations) -> %.1f s total\n', ...
        id.time/max(id.nIter,1), id.nIter, id.time);
fprintf('  ratio, both implementations optimised : %.1fx\n', ig.time/max(id.time,eps));
fprintf('  ratio, literal Alg.1 vs GD            : %.0fx   (paper: ~250x)\n', ...
        naiveIter*max(ig.nIter,1)/max(id.time,eps));
fprintf('  -> the 250x is a property of the implementation (cost per candidate),\n');
fprintf('     not of the algorithms: same Algorithm 1, same iterates.\n');
fprintf('  wall-clock on this machine: greedy %.1f s, GD %.1f s (paper: 1.05 s for GD)\n', ...
        ig.time, id.time);

%% ---------------- Part 4 : where does "+2.26 dB" come from ------------
fprintf('\n--- Part 4: the stated +2.26/+2.11 dB needs a baseline of ~20.4 dB ---\n');
fprintf('  %-11s |  u1     u2     u3     u4   |  mean  | greedy | delta\n', 'baseline');
conv = struct('name',{},'ls',{},'gr',{},'delta',{});
for nowrap = [false true]
    q = p;  q.nowrap = nowrap;
    blq = jp_block(q, alpha);
    [ph, ta] = jp_design(q, az, el, alpha, 'joint-ls');
    R = jp_eval(q, ph, ta, az, el, blq);
    o = oG; o.mode = 'joint';
    [pg, tg] = jp_greedy(q, az, el, alpha, o);
    Rg = jp_eval(q, pg, tg, az, el, blq);
    if nowrap, nm = 'raw nu'; else, nm = 'unwrapped'; end
    fprintf('  %-11s | %6.3f %6.3f %6.3f %6.3f | %6.3f | %6.3f | %+5.2f\n', ...
            nm, R.gdb, mean(R.gdb), mean(Rg.gdb), mean(Rg.gdb)-mean(R.gdb));
    conv(end+1) = struct('name',nm,'ls',mean(R.gdb),'gr',mean(Rg.gdb), ...
                         'delta',mean(Rg.gdb)-mean(R.gdb));  %#ok<SAGROW>
end
fprintf('  needed for +2.26 dB : %.2f dB  (= joint greedy %.2f - 2.26)\n', ...
        mean(G5a(2,:))-2.26, mean(G5a(2,:)));
fprintf('  my separated LS is %.2f dB, i.e. joint greedy - sep LS = %+.2f dB\n', ...
        mean(G5a(4,:)), mean(G5a(2,:))-mean(G5a(4,:)));
fprintf('  (Fig.2a check: at alpha=0.10 the unwrapped LS gives 21.97/18.63 vs paper 22.1/18.1;\n');
fprintf('   the raw-nu convention gives -3.05/-4.60, so the paper does use an unwrapped target.)\n');

%% ---------------- save ----------------
save(fullfile(out,'a6_greedy_gd.mat'),'names','G5a','T5a','ITS','avec','lbl2','G5b','T5b','ITb', ...
     'chk','ig','id','igs','ids','tE','naiveIter','ZETA_G','pap5a','conv');

fid = fopen(fullfile(out,'a6_fig5a.csv'),'w');
fprintf(fid,'design,u1,u2,u3,u4,mean,iters,time_s\n');
for k = 1:numel(names)
    fprintf(fid,'%s,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.3f\n', names{k}, G5a(k,:), mean(G5a(k,:)), ITS(k), T5a(k));
end
fclose(fid);

fid = fopen(fullfile(out,'a6_fig5b.csv'),'w');
fprintf(fid,'alpha,LSjoint,GrJoint,GDJoint,LSsep,GrSep,GDSep,tGr_s,tGD_s\n');
for k = 1:numel(avec)
    fprintf(fid,'%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f\n', avec(k), G5b(:,k), T5b(2,k), T5b(3,k));
end
fclose(fid);

fid = fopen(fullfile(out,'a6_summary.txt'),'w');
fprintf(fid,'A6 (Fig.5): greedy Alg.1 + Adam GD   [nanyangnorth]\n\n');
fprintf(fid,'equivalence: incremental vs literal Alg.1, joint diff %.2e, sep diff %.2e (Gl units)\n', ...
        chk.joint.diff, chk.sep.diff);
fprintf(fid,'Fig5a per-user max |mine-paper|: greedy j %.2f, GD j %.2f, greedy s %.2f, GD s %.2f dB\n', ...
        max(abs(G5a(2,:)-pap5a.jointGreedy)), max(abs(G5a(3,:)-pap5a.jointGD)), ...
        max(abs(G5a(5,:)-pap5a.sepGreedy)), max(abs(G5a(6,:)-pap5a.sepGD)));
fprintf(fid,'Fig5a joint delta over joint-LS: greedy %+0.2f, GD %+0.2f (paper +2.26 / +2.11)\n', ...
        mean(G5a(2,:))-mean(G5a(1,:)), mean(G5a(3,:))-mean(G5a(1,:)));
fprintf(fid,'Fig5b mean|GD-greedy| joint %.3f dB, sep %.3f dB (paper: comparable)\n', ...
        mean(abs(G5b(3,:)-G5b(2,:))), mean(abs(G5b(6,:)-G5b(5,:))));
fprintf(fid,'runtime: greedy %.2f s (%.2f s/pass, %d passes), GD %.2f s (%.3f s/it, %d it)\n', ...
        ig.time, ig.time/ig.nIter, ig.nIter, id.time, id.time/id.nIter, id.nIter);
fprintf(fid,'ratio optimised %.1fx; literal Alg.1 vs GD %.0fx (paper ~250x)\n', ...
        ig.time/id.time, naiveIter*ig.nIter/id.time);
fprintf(fid,'baseline needed for +2.26 dB: %.2f dB; my joint LS %.2f, sep LS %.2f\n', ...
        mean(G5a(2,:))-2.26, mean(G5a(1,:)), mean(G5a(4,:)));
fclose(fid);

fprintf('\nwrote a6_* into %s\n', out);
end
