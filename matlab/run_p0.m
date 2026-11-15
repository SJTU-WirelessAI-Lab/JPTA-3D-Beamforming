function run_p0()
% P0 self-check of this paper's core.  All assertions must pass before any
% experiment is allowed to run (frozen-kernel policy of the repro plan).

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));

fprintf('===== P0 SELF-CHECK : this paper core =====\n');
fprintf('matlab %s on %s\n', version, computer);
p = jp_params();
fprintf('Naz x Nel = %d x %d, Nelem = %d, Gmax = %.3f dB\n\n', p.Naz, p.Nel, p.Nelem, p.Gmax_dB);

m = jp_mvec(p);  f = p.fc + m*p.df;
pass = true;  T = {};

%% --- check 1: Nu = 1, alpha = 1  ->  delays must be constant (pure PAA) ---
az = -60; el = 90;                       % a sector-edge direction (non-trivial phases)
[phi, tau, info] = jp_design(p, az, el, 1, 'joint-ls');
g = 10*log10(jp_gain_freqs(p, phi, tau, az, el, f));
c1 = (max(abs(tau(:))) < 1e-9) && (abs(mean(g)-p.Gmax_dB) < 0.2) && (max(g)-min(g) < 0.5);
fprintf('[1] single user  : max|tau| = %.2e s | mean gain = %.3f dB | ripple = %.3f dB\n', ...
        max(abs(tau(:))), mean(g), max(g)-min(g));
fprintf('    required delay span = %.2f ns\n', info.tau_span_raw*1e9);
fprintf('    -> %s\n', tf(c1));
pass = pass && c1;
T(end+1,:) = {'P0-1 single-user degenerates to PAA', c1, sprintf('tau=%.1e, gain=%.3f dB, ripple=%.3f dB', max(abs(tau(:))), mean(g), max(g)-min(g))};

%% --- check 2: no solution may exceed 10*log10(Naz*Nel) ---
worst = -inf;  worstcase = '';
for Nu = [1 2 4 5]
    [azv, elv] = jp_users(Nu);
    al = ones(1,Nu)/Nu;
    modes = {'joint-ls','sep-ls','joint-inf','sep-inf'};
    for mi = 1:numel(modes)
        [phi, tau] = jp_design(p, azv, elv, al, modes{mi});
        R = jp_eval(p, phi, tau, azv, elv, jp_block(p, al));
        if max(R.gdb) > worst
            worst = max(R.gdb);  worstcase = sprintf('Nu=%d %s', Nu, modes{mi});
        end
    end
end
c2 = worst <= p.Gmax_dB + 1e-6;
fprintf('[2] gain bound    : max observed = %.4f dB (limit %.4f dB) at %s\n', worst, p.Gmax_dB, worstcase);
fprintf('    -> %s\n', tf(c2));
pass = pass && c2;
T(end+1,:) = {'P0-2 gain <= 10log10(Naz*Nel)', c2, sprintf('max %.4f dB vs limit %.4f dB', worst, p.Gmax_dB)};

%% --- check 3: with Nel = 1 the separated design must equal the joint one ---
q = p;  q.Nel = 1;  q.Nelem = q.Naz*q.Nel;  q.Gmax_dB = 10*log10(q.Nelem);
azv = [-60 -20 20 55]; elv = 90*ones(1,4);        % elevation is degenerate when Nel=1
al  = ones(1,4)/4;
[phiA, tauA] = jp_design(q, azv, elv, al, 'joint-ls');
[phiB, tauB] = jp_design(q, azv, elv, al, 'sep-ls');
RA = jp_eval(q, phiA, tauA, azv, elv, jp_block(q, al));
RB = jp_eval(q, phiB, tauB, azv, elv, jp_block(q, al));
d3 = max(abs(RA.gdb - RB.gdb));
c3 = (d3 < 1e-9) && (max(abs(tauA(:)-tauB(:))) < 1e-12);
fprintf('[3] Nel=1 ULA    : max|joint-separated| = %.2e dB, max|dtau| = %.2e s\n', d3, max(abs(tauA(:)-tauB(:))));
fprintf('    -> %s\n', tf(c3));
pass = pass && c3;
T(end+1,:) = {'P0-3 ULA: joint == separated', c3, sprintf('dGain=%.2e dB', d3)};

%% --- check 4: delay fits the hardware range for the default scenarios ---
ok = true;  spans = [];
for Nu = [2 4 8 16]
    [azv, elv] = jp_users(Nu);
    [~, ~, inf4] = jp_design(p, azv, elv, ones(1,Nu)/Nu, 'joint-ls');
    spans(end+1) = inf4.tau_span_raw*1e9; %#ok<AGROW>
    ok = ok && (inf4.tau_span_raw <= p.tau_max);
end
c4 = ok;
fprintf('[4] delay range  : required spans (Nu=2,4,8,16) = %s ns, tau_max = %.0f ns\n', ...
        mat2str(round(spans,1)), p.tau_max*1e9);
fprintf('    -> %s\n', tf(c4));
pass = pass && c4;
T(end+1,:) = {'P0-4 delay within tau_max', c4, sprintf('spans %s ns', mat2str(round(spans,1)))};

%% --- summary ---
fprintf('\n===== SUMMARY: %s =====\n', upper(tf(pass)));
for i = 1:size(T,1)
    fprintf('  %-36s : %-4s  (%s)\n', T{i,1}, tf(T{i,2}), T{i,3});
end

out = fullfile(here,'out'); if ~exist(out,'dir'), mkdir(out); end
fid = fopen(fullfile(out,'p0_selfcheck.txt'),'w');
for i = 1:size(T,1)
    fprintf(fid, '%s\t%d\t%s\n', T{i,1}, T{i,2}, T{i,3});
end
fprintf(fid, 'ALL_PASS\t%d\n', pass);
fclose(fid);
fprintf('\nwrote %s\n', fullfile(out,'p0_selfcheck.txt'));
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
