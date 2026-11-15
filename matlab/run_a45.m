function run_a45()
% A4 / A5 of the repro plan = C16 Figs. 3 and 4.
%   A4: max-over-subcarriers gain in the (azimuth, elevation) plane for the
%       separated and the joint LS design, Nu = 5, equal BW allocation.
%   A5: gain in the (azimuth, frequency) plane at theta_el = 105 deg for the
%       joint LS design, with user 3's band highlighted.
%
% Cross-check built into A5: with Nu=5 and alpha=0.2 each, user 3 (i=2) owns
% m = -79..79, i.e. +-9.48 MHz - the paper's Fig. 4 says [-9.5, 9.5) MHz.
% This validates the subcarrier-allocation convention of jp_block.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'core'));
addpath(fullfile(here,'alg'));
out = fullfile(here,'out'); if ~exist(out,'dir'), mkdir(out); end

p  = jp_params_c16();
Nu = 5;
[az, el] = jp_users_c16(Nu);
a  = ones(1,Nu)/Nu;
bl = jp_block(p, a);
mm = jp_mvec(p);  fm = p.fc + mm*p.df;

i3 = 3;
sel3 = (bl == i3);
fprintf('===== A4 / A5 : C16 Figs. 3 and 4 (Nu=%d, equal BW) =====\n', Nu);
fprintf('user 3 band: m = %d..%d  ->  [%.2f, %.2f] MHz\n', ...
        mm(find(sel3,1,'first')), mm(find(sel3,1,'last')), ...
        mm(find(sel3,1,'first'))*p.df/1e6, mm(find(sel3,1,'last'))*p.df/1e6);

% ---------------- A5 : azimuth vs frequency at el = 105 deg ----------------
[phiJ, tauJ] = jp_design(p, az, el, a, 'joint-ls');
azv = -60:0.25:60;
H = zeros(numel(azv), numel(fm));
for k = 1:numel(azv)
    H(k,:) = 10*log10(jp_gain_freqs(p, phiJ, tauJ, azv(k), 105, fm));
end

fg = figure('Visible','off','Position',[100 100 760 460]);
foff = (fm - p.fc)/1e6;
imagesc(foff, azv, H); axis xy; hold on;
cb = colorbar; cb.Label.String = 'Beamforming gain (dB)';
plot([-9.48 -9.48], [-60 60], 'r-', 'LineWidth',2);
plot([ 9.48  9.48], [-60 60], 'r-', 'LineWidth',2);
xlabel('Subcarrier frequency offset (MHz)');
ylabel('Azimuth (deg)');
title(sprintf('C16 Fig.4 reproduction: joint LS, \\theta_{el}=105^\\circ, user 3 at (0^\\circ,105^\\circ)'));
caxis([8 26]); xlim([-48 48]);
print(fg, fullfile(out,'a5_az_freq_joint.png'), '-dpng', '-r180');
close(fg);

% quantitative check at the user's own direction (az = 0)
k0 = find(abs(azv) < 1e-9, 1);
fprintf('A5: at az = 0 deg, mean gain in band = %.2f dB, outside band = %.2f dB (Gmax %.2f)\n', ...
        mean(H(k0,sel3)), mean(H(k0,~sel3)), p.Gmax_dB);

% ---------------- A4 : azimuth vs elevation, max over subcarriers ----------
fb = fm(1:8:end);                       % decimated grid is enough for a max-hold map
azg = -60:1:60;  elg = 90:0.5:120;
modes = {'sep-ls','joint-ls'};
tag   = {'sep','joint'};
for mi = 1:2
    [ph, ta] = jp_design(p, az, el, a, modes{mi});
    Mx = zeros(numel(azg), numel(elg));
    for e = 1:numel(elg)
        for k = 1:numel(azg)
            g = jp_gain_freqs(p, ph, ta, azg(k), elg(e), fb);
            Mx(k,e) = 10*log10(max(g));
        end
    end
    fg = figure('Visible','off','Position',[100 100 700 470]);
    imagesc(elg, azg, Mx); axis xy; hold on;
    cb = colorbar; cb.Label.String = 'max gain over subcarriers (dB)';
    plot(el, az, 'rd', 'MarkerSize',10, 'MarkerFaceColor','r');
    xlabel('Elevation (deg)'); ylabel('Azimuth (deg)');
    title(sprintf('C16 Fig.3 reproduction (%s LS), N_u=%d, red = desired users', tag{mi}, Nu));
    caxis([5 26]);
    print(fg, fullfile(out, sprintf('a4_pattern_%s.png', tag{mi})), '-dpng', '-r180');
    close(fg);
    save(fullfile(out, sprintf('a4_%s.mat', tag{mi})), 'Mx','azg','elg');
    fprintf('A4(%s): done, size %dx%d\n', tag{mi}, size(Mx,1), size(Mx,2));
end
fprintf('\nfigures written to %s\n', out);
end
