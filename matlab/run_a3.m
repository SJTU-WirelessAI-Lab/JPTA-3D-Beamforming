function run_a3()
% A3 of the repro plan = Fig. 2(b): per-user average gain for Nu = 5 with
% an unequal bandwidth allocation alpha = [0.3 0.2 0.15 0.1 0.25].
% Paper's claims (max gain difference among users):
%   joint inf-norm 3.0 dB   vs  joint LS 5.0 dB
%   sep   inf-norm 3.2 dB   vs  sep   LS 5.4 dB

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));
out = fullfile(here,'out'); if ~exist(out,'dir'), mkdir(out); end

p  = jp_params();
Nu = 5;
[az, el] = jp_users(Nu);
a  = [0.3 0.2 0.15 0.1 0.25];
bl = jp_block(p, a);

modes = {'joint-inf','joint-ls','sep-inf','sep-ls'};
lbl   = {'Joint Inf Norm','Joint LS [10]','Separated Inf Norm','Separated LS [10]'};
claim = [3.0 5.0 3.2 5.4];

fprintf('===== A3 : Fig.2(b), Nu = %d, alpha = [%s] =====\n', Nu, num2str(a));
fprintf('users: az = [%s], el = [%s] deg\n\n', num2str(az), num2str(el));

gdb = zeros(numel(modes), Nu);
for mi = 1:numel(modes)
    [phi, tau] = jp_design(p, az, el, a, modes{mi});
    R = jp_eval(p, phi, tau, az, el, bl);
    gdb(mi,:) = R.gdb;
end

fprintf('  mode            |   U1     U2     U3     U4     U5   | spread | paper\n');
fprintf('  ----------------+------------------------------------+--------+------\n');
for mi = 1:numel(modes)
    sp = max(gdb(mi,:)) - min(gdb(mi,:));
    fprintf('  %-15s | %s | %5.2f  | %4.1f\n', lbl{mi}, ...
            sprintf('%6.2f ', gdb(mi,:)), sp, claim(mi));
end

% ---- figure (Fig. 2b) ----
fg = figure('Visible','off','Position',[100 100 640 480]);
hold on; box on; grid on;
st = {'-o','-s','--^','--d'};
cl = [0.85 0.33 0.10; 0 0.45 0.74; 0.47 0.67 0.19; 0.49 0.18 0.56];
for mi = 1:numel(modes)
    plot(1:Nu, gdb(mi,:), st{mi}, 'Color', cl(mi,:), 'LineWidth',1.6, ...
         'MarkerSize',7, 'MarkerFaceColor', cl(mi,:));
end
xlabel('User index'); ylabel('Average gain (dB)');
title(sprintf('Fig.2(b) reproduction: N_u=%d, \\alpha=[%s]', Nu, num2str(a)));
legend(lbl,'Location','northeast');
xlim([0.8 5.2]); ylim([15 26]);
print(fg, fullfile(out,'a3_fairness.png'), '-dpng', '-r180');
close(fg);
save(fullfile(out,'a3_fairness.mat'),'gdb','az','el','a','lbl');
fprintf('\nwrote %s\n', fullfile(out,'a3_fairness.png'));
end
