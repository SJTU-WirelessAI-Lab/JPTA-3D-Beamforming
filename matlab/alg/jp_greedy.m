function [phi, tau, info] = jp_greedy(p, az, el, alpha, opt)
% Greedy beam design , Sec. III-C, i.e. Algorithm 1 (joint) and its
% separated counterpart.
%
%   opt.mode     'joint' (default) | 'sep'
%   opt.zeta     convergence threshold, |Gl_new - Gl_old| < zeta*|Gl_new|
%   opt.maxIter  hard iteration cap
%
% Search spaces are exactly those of Algorithm 1:
%   tau grid  = 0 : tau_p : tau_max        (line 4)
%   phi grid  = 0 : 2*pi/2^beta : 2*pi     (line 5)
% and the objective is Gl of eq. (5).
%
% IMPLEMENTATION NOTE (this is the only deviation from a literal reading):
% updating a single element (y,z) only changes that element's term in the
% array sum, so the candidate objective is evaluated from
%       S_rest = S - contribution_of_(y,z)
% instead of recomputing the whole gain for every candidate.  The result is
% numerically identical to the naive implementation (checked against
% jp_greedy_naive.m); only the cost changes, from O(N_elem) per candidate to
% O(1) per candidate.  The runtime claim of the paper (GD ~250x faster than
% greedy) therefore depends on this implementation choice, not only on the
% algorithms - see run_a6.m.
%
% Joint   : one axis of N_elem elements, external weight Wext = 1.
% Separated: two axes (Naz, then Nel); the other axis' |S|^2 is Wext.

if nargin < 5, opt = struct(); end
opt = jp_filldef(opt, struct('mode','joint','zeta',1e-4,'maxIter',40,'verbose',false));

blk = jp_block(p, alpha);
m   = jp_mvec(p);   f = p.fc + m*p.df;   Nk = p.M+1;
Nu  = numel(az);

selu = false(Nk, Nu);  ns = zeros(1,Nu);
for i = 1:Nu
    selu(:,i) = (blk(:) == i);
    ns(i) = nnz(selu(:,i));
end
if any(ns == 0)
    error('jp_greedy: a user was allocated 0 subcarriers');
end

tg = (0:p.tau_p:p.tau_max).';            % delay candidates
q  = 2*pi/2^p.beta_ps;
pg = (0:q:(2*pi-q)).';                   % phase candidates
D  = exp(1j*2*pi*(tg*f));                % delay phase factor, nTg x Nk

cfg = struct('f',f,'D',D,'tg',tg,'pg',pg,'selu',selu,'ns',ns,'N',p.Nelem);

info.mode = opt.mode;  info.zeta = opt.zeta;

switch lower(opt.mode)

    % ================= joint =========================================
    case 'joint'
        Om = zeros(p.Nelem, Nk, Nu);
        for i = 1:Nu
            nu = jp_nu(p, az(i), el(i));
            Om(:,:,i) = nu(:) * (f/p.fc);
        end
        Mx = exp(-1j*Om);                        % steering factor

        [phi, tau] = jp_design(p, az, el, alpha, 'joint-ls');
        phix = mod(round(phi(:)/q)*q, 2*pi);
        taux = min(max(round(tau(:)/p.tau_p)*p.tau_p, 0), p.tau_max);
        idxT = round(taux/p.tau_p) + 1;

        Hx = phix + 2*pi*(taux*f);               % N_elem x Nk
        A  = exp(1j*Hx) .* Mx;                   % N_elem x Nk x Nu
        Sx = reshape(sum(A,1), Nk, Nu);         % Nk x Nu
        Wext = ones(Nk, Nu);

        Gl = globj(Sx, Wext, selu, ns, p.Nelem);
        Glhist = Gl;
        t0 = tic;

        for it = 1:opt.maxIter
            Gl0 = Gl;
            [Sx, phix, taux, idxT, Hx, A, ~] = sweep_axis(A, Sx, phix, taux, ...
                idxT, Hx, Mx, cfg, Wext);
            Gl1 = globj(Sx, Wext, selu, ns, p.Nelem);
            Glhist(end+1) = Gl1;                 %#ok<AGROW>
            if abs(Gl1 - Gl0) < opt.zeta*abs(Gl1)
                Gl = Gl1;  break;
            end
            Gl = Gl1;
            if opt.verbose
                fprintf('    it %2d : Gl = %.4f  (d = %+.4f)\n', it, Gl, Gl-Gl0);
            end
        end
        info.nIter = numel(Glhist) - 1;
        info.time  = toc(t0);
        info.Gl    = Gl;
        info.Glhist = Glhist;

        phi = mod(reshape(phix, p.Naz, p.Nel), 2*pi);
        tau = reshape(taux, p.Naz, p.Nel);

    % ================= separated =====================================
    case 'sep'
        [~, ~, s] = jp_design(p, az, el, alpha, 'sep-ls');
        phia = mod(round(s.phi_az(:)/q)*q, 2*pi);
        taua = min(max(round(s.tau_az(:)/p.tau_p)*p.tau_p, 0), p.tau_max);
        phie = mod(round(s.phi_el(:)/q)*q, 2*pi);
        taue = min(max(round(s.tau_el(:)/p.tau_p)*p.tau_p, 0), p.tau_max);

        OmA = zeros(p.Naz, Nk, Nu);  OmE = zeros(p.Nel, Nk, Nu);
        for i = 1:Nu
            a = az(i)*pi/180;  e = el(i)*pi/180;
            OmA(:,:,i) = (pi*(0:p.Naz-1).'*(sin(a)*sin(e))) * (f/p.fc);
            OmE(:,:,i) = (pi*(0:p.Nel-1).'*cos(e))            * (f/p.fc);
        end
        MA = exp(-1j*OmA);  ME = exp(-1j*OmE);

        HA = phia + 2*pi*(taua*f);
        HE = phie + 2*pi*(taue*f);
        AA = exp(1j*HA) .* MA;    Saz = reshape(sum(AA,1), Nk, Nu);
        AE = exp(1j*HE) .* ME;    Sel = reshape(sum(AE,1), Nk, Nu);

        idxTa = round(taua/p.tau_p) + 1;
        idxTe = round(taue/p.tau_p) + 1;

        Gl = globj(Saz, abs(Sel).^2, selu, ns, p.Nelem);
        Glhist = Gl;
        t0 = tic;

        for it = 1:opt.maxIter
            Gl0 = Gl;
            % ---- azimuth axis (elevation gains held fixed) ----
            cfg.tauCap = p.tau_max - max(taue);
            [Saz, phia, taua, idxTa, HA, AA, ~] = sweep_axis(AA, Saz, phia, taua, ...
                idxTa, HA, MA, cfg, abs(Sel).^2);

            % ---- elevation axis (azimuth gains held fixed) ----
            cfg.tauCap = p.tau_max - max(taua);
            [Sel, phie, taue, idxTe, HE, AE, ~] = sweep_axis(AE, Sel, phie, taue, ...
                idxTe, HE, ME, cfg, abs(Saz).^2);

            Gl1 = globj(Saz, abs(Sel).^2, selu, ns, p.Nelem);
            Glhist(end+1) = Gl1;                 %#ok<AGROW>
            if abs(Gl1 - Gl0) < opt.zeta*abs(Gl1)
                Gl = Gl1;  break;
            end
            Gl = Gl1;
            if opt.verbose
                fprintf('    it %2d : Gl = %.4f  (d = %+.4f)\n', it, Gl, Gl-Gl0);
            end
        end
        info.nIter = numel(Glhist) - 1;
        info.time  = toc(t0);
        info.Gl    = Gl;
        info.Glhist = Glhist;

        phi = mod(phia + phie.', 2*pi);
        tau = taua + taue.';
        if max(tau(:)) > p.tau_max + eps
            tau = min(tau, p.tau_max);           % same clipping as jp_design
        end

    otherwise
        error('jp_greedy: unknown mode %s', opt.mode);
end

info.tau_span_raw = max(tau(:));
info.saturated    = max(tau(:)) > p.tau_max + eps;
end

% =====================================================================
function [Sx, phix, taux, idxT, Hx, A, Gl] = sweep_axis(A, Sx, phix, taux, ...
        idxT, Hx, Mx, cfg, Wext)
% One full pass of Algorithm 1 over ALL elements of a single axis: first the
% delay grid for every element (lines 11-13), then the phase grid for every
% element (lines 14-16).
Nk = size(Sx,1);  Nu = size(Sx,2);
N  = cfg.N;  f = cfg.f;  D = cfg.D;  tg = cfg.tg;  pg = cfg.pg;
nT = numel(tg);  nP = numel(pg);
hasCap = isfield(cfg,'tauCap');

for e = 1:numel(phix)
    Oe   = reshape(Mx(e,:,:), Nk, Nu);
    rest = Sx - reshape(A(e,:,:), Nk, Nu);

    % ---- delay sweep for element e ----
    Uc  = exp(1j*phix(e)) .* D;                    % nT x Nk
    obj = cand_obj(Uc, rest, Wext, Mx, e, cfg, Nk, Nu);
    if hasCap, obj(tg > cfg.tauCap + 1e-15) = -Inf; end
    [~, kb] = max(obj);
    taux(e) = tg(kb);  idxT(e) = kb;
    Hx(e,:) = phix(e) + 2*pi*taux(e)*f;

    % ---- phase sweep for element e ----
    Uc  = exp(1j*pg) .* D(idxT(e),:);              % nP x Nk
    obj = cand_obj(Uc, rest, Wext, Mx, e, cfg, Nk, Nu);
    [~, kb] = max(obj);
    phix(e) = pg(kb);
    Hx(e,:) = phix(e) + 2*pi*taux(e)*f;

    Ae = exp(1j*Hx(e,:)) .* Mx(e,:,:);             % 1 x Nk x Nu
    A(e,:,:) = Ae;
    Sx = rest + reshape(Ae, Nk, Nu);
end
Gl = globj(Sx, Wext, cfg.selu, cfg.ns, N);
end

% ---------------------------------------------------------------------
function obj = cand_obj(Uc, rest, Wext, Mx, e, cfg, Nk, Nu)
% Gl of eq. (5) for every candidate row of Uc (= exp(j*h_cand)), with the
% other elements held fixed in rest.
obj = zeros(size(Uc,1), 1);
N   = cfg.N;
for i = 1:Nu
    Cc = Uc .* reshape(Mx(e,:,i), 1, Nk);          % nCand x Nk
    Sn = Cc + rest(:,i).';                         % nCand x Nk
    obj = obj + 10*log10( sum((abs(Sn).^2).*Wext(:,i).'.*cfg.selu(:,i).', 2) ...
                          / (N*cfg.ns(i)) );
end
end

% ---------------------------------------------------------------------
function Gl = globj(Sx, Wext, selu, ns, N)
% Gl = sum_i 10*log10( mean_{k in F_i} |Sx(k,i)|^2 * Wext(k,i) / N )
Nu = size(Sx,2);
Gl = 0;
for i = 1:Nu
    Gm = sum((abs(Sx(:,i)).^2).*Wext(:,i).*selu(:,i)) / (N*ns(i));
    if ~(Gm > 0), Gm = realmin; end
    Gl = Gl + 10*log10(Gm);
end
end
