function run_a6_conv()
% Diagnostic for A6: (A) does the paper's "+2.26 dB" correspond to a
% different phase-unwrapping convention in the analytical starting point?
% (B) smoke test of the bounded-parameterisation GD.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));
p = jp_params_c16();

fprintf('=== A: LS convention, joint, Nu = 4 equal BW ===\n');
fprintf('  convention |  u1     u2     u3     u4   |  mean  | greedy mean | delta\n');
for nowrap = [false true]
    q = p;  q.nowrap = nowrap;
    Nu = 4;  [az, el] = jp_users_c16(Nu);  al = ones(1,Nu)/Nu;
    blk = jp_block(q, al);
    [ph, ta] = jp_design(q, az, el, al, 'joint-ls');
    R = jp_eval(q, ph, ta, az, el, blk);
    o = struct('mode','joint','zeta',1e-4,'maxIter',40);
    [pg, tg, ig] = jp_greedy(q, az, el, al, o);
    Rg = jp_eval(q, pg, tg, az, el, blk);
    if nowrap, nm = 'raw nu    '; else, nm = 'unwrapped '; end
    fprintf('  %s | %6.3f %6.3f %6.3f %6.3f | %6.3f | %6.3f      | %+5.2f (greedy %+5.2f)\n', ...
            nm, R.gdb, mean(R.gdb), mean(Rg.gdb), mean(Rg.gdb)-mean(R.gdb), ...
            mean(Rg.gdb)-mean(R.gdb));
    fprintf('    tau span LS %.1f ns, greedy %.1f ns, greedy iters %d\n', ...
            max(ta(:))*1e9, max(tg(:))*1e9, ig.nIter);
end

fprintf('\n=== A2: LS convention, Nu = 2, alpha = 0.10 (paper Fig.2a: 22.1 / 18.1) ===\n');
Nu = 2;  [az2, el2] = jp_users_c16(Nu);  a2 = [0.10 0.90];
blk2 = jp_block(p, a2);
for nowrap = [false true]
    q = p;  q.nowrap = nowrap;
    [ph, ta] = jp_design(q, az2, el2, a2, 'joint-ls');
    R = jp_eval(q, ph, ta, az2, el2, blk2);
    [ph, ta] = jp_design(q, az2, el2, a2, 'sep-ls');
    Rs = jp_eval(q, ph, ta, az2, el2, blk2);
    if nowrap, nm = 'raw nu   '; else, nm = 'unwrapped'; end
    fprintf('  %s : joint LS %.2f (paper 22.1) | sep LS %.2f (paper 18.1)\n', ...
            nm, R.Gl_avg, Rs.Gl_avg);
end

fprintf('\n=== B: GD smoke test, Nu = 4 equal BW (greedy target 22.709) ===\n');
Nu = 4;  [az, el] = jp_users_c16(Nu);  al = ones(1,Nu)/Nu;
blk = jp_block(p, al);
[ph, ta] = jp_design(p, az, el, al, 'joint-ls');
R = jp_eval(p, ph, ta, az, el, blk);
fprintf('  LS start: mean %.3f dB, span %.1f ns\n', mean(R.gdb), max(ta(:))*1e9);
for mxi = [100 300 600]
    [pg, tg, id] = jp_gd(p, az, el, al, struct('mode','joint','maxIter',mxi));
    fprintf('  GD joint (cap %3d): mean %.3f dB, iters %3d, %5.2f s, span %5.1f ns, Gl_drop %+.1e\n', ...
            mxi, mean(id.gdb), id.nIter, id.time, id.tau_span_raw*1e9, id.Gl_drop);
end
for mxi = [300 600]
    [pg, tg, ids] = jp_gd(p, az, el, al, struct('mode','sep','maxIter',mxi));
    fprintf('  GD sep   (cap %3d): mean %.3f dB, iters %3d, %5.2f s, span %5.1f ns, Gl_drop %+.1e\n', ...
            mxi, mean(ids.gdb), ids.nIter, ids.time, ids.tau_span_raw*1e9, ids.Gl_drop);
end
end
