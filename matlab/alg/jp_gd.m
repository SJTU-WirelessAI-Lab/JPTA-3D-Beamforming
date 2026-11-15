function [phi, tau, info] = jp_gd(p, az, el, alpha, opt)
% Gradient-descent beam design , Sec. III-D (Adam, lr = 0.1), for the
% joint and the separated parameterisation.
%
%   opt.mode      'joint' (default) | 'sep'
%   opt.lr        0.1      (paper: chosen by grid search)
%   opt.zeta      1e-4     (same convergence test as Algorithm 1, line 18)
%   opt.maxIter   600
%   opt.quantize  false
%
% ---------------------------------------------------------------------------
% PARAMETERISATION (documented deviations, gaps G11/G12 of the repro plan)
%
% (i) The paper differentiates w.r.t. (phi, tau) with the SAME learning rate
%     0.1.  Under Adam the step in a parameter is ~lr regardless of gradient
%     scale, so lr = 0.1 on tau *in seconds* would move the per-subcarrier
%     phase slope b = 2*pi*df*tau by 2*pi*120e3*0.1 ~ 7.5e4 rad per step,
%     while lr = 0.1 on tau *in units of tau_p = 2.5 ns* moves it by 1.9e-3
%     rad per step (800 steps to cross the range).  The only scale on which a
%     single lr = 0.1 is meaningful for both variables is
%         a = phi + 2*pi*fc*tau   (phase at the centre subcarrier)   [rad]
%         b = 2*pi*df*tau         (phase slope per subcarrier)       [rad]
%     i.e. lr = 0.1 rad.  We optimise in (a,b).
%
% (ii) Hardware constraint: tau in [0, tau_max] implies b in [0, Bmax] with
%     Bmax = 2*pi*df*tau_max = 0.151 rad, a window only 2.4% of the 2*pi
%     period of b.  An unconstrained Adam on b leaves that window within a
%     few steps - verified: it ends pinned at the tau_max clip and the
%     resulting design collapses once mapped back (Gl drops by ~70 dB).
%     We therefore optimise a bounded variable
%         b = Bmax * sigmoid(u)
%     so every iterate is feasible and the mapping back is exact
%     (info.Gl_drop ~ 0).  A common shift of all b is gain-invariant (it
%     multiplies the array sum by exp(j*delta*m)), so pinning the lower edge
%     at 0 costs nothing.  For the separated case tau = tau_az + tau_el, so
%     each axis gets half of the range.
% ---------------------------------------------------------------------------

if nargin < 5, opt = struct(); end
opt = jp_filldef(opt, struct('mode','joint','lr',0.1,'zeta',1e-5,'win',20, ...
                             'maxIter',600,'quantize',false,'verbose',false));

blk = jp_block(p, alpha);
m   = jp_mvec(p);   f = p.fc + m*p.df;   Nk = p.M+1;
Nu  = numel(az);    N  = p.Nelem;
mrow = m(:).';

selu = false(Nk, Nu);  ns = zeros(1,Nu);
for i = 1:Nu
    selu(:,i) = (blk(:) == i);
    ns(i) = nnz(selu(:,i));
end

Bmax = 2*pi*p.df*p.tau_max;              % 0.151 rad
b1 = 0.9;  b2 = 0.999;  epsA = 1e-8;     % Adam hyper-parameters
info.mode = opt.mode;
info.Bmax = Bmax;

lg = @(z) log(z./(1-z));
sg = @(z) 1./(1+exp(-z));

switch lower(opt.mode)

    % ================= joint =========================================
    case 'joint'
        [phi0, tau0] = jp_design(p, az, el, alpha, 'joint-ls');
        a  = phi0(:) + 2*pi*p.fc*tau0(:);
        b0 = 2*pi*p.df*tau0(:);
        u  = lg(min(max(b0/Bmax, 1e-3), 1-1e-3));

        Om = zeros(N, Nk, Nu);
        for i = 1:Nu
            nu = jp_nu(p, az(i), el(i));
            Om(:,:,i) = nu(:) * (f/p.fc);
        end
        Mx = exp(-1j*Om);

        ma = zeros(N,1); va = zeros(N,1);
        mu = zeros(N,1); vu = zeros(N,1);
        tt = 0;  Glhist = [];  b = Bmax*sg(u);
        t0 = tic;

        for it = 1:opt.maxIter
            s  = sg(u);  b = Bmax*s;
            H  = a + b*mrow;
            A  = exp(1j*H) .* Mx;
            S  = reshape(sum(A,1), Nk, Nu);
            Gk = abs(S).^2 / N;
            Gmean = sum(Gk.*selu, 1)./ns;
            Gl = sum(10*log10(Gmean));
            Glhist(end+1) = Gl;                    %#ok<AGROW>
            if converged(Glhist, opt.zeta, opt.win), break; end

            w = 10/log(10) ./ Gmean;
            Q = zeros(Nk, Nu);
            for i = 1:Nu
                Q(selu(:,i),i) = w(i)*(2/N)/ns(i)*conj(S(selu(:,i),i));
            end
            W = zeros(N, Nk);
            for i = 1:Nu
                W = W + A(:,:,i).*Q(:,i).';
            end
            ga = -imag(sum(W, 2));
            gb = -imag(sum(W.*mrow, 2));

            [a, ma, va, tt] = adam(a, ga, ma, va, tt, opt.lr, b1, b2, epsA);
            % chain rule through b = Bmax*sigmoid(u)
            [u, mu, vu, ~ ] = adam(u, gb.*(Bmax*s.*(1-s)), mu, vu, tt, opt.lr, b1, b2, epsA);

            if opt.verbose && mod(it,50) == 0
                fprintf('    gd it %3d : Gl = %.4f  tau span %.1f ns\n', it, Gl, ...
                        (max(b)-min(b))/(2*pi*p.df)*1e9);
            end
        end
        info.nIter = numel(Glhist)-1;
        info.Gl_opt = Gl;

        s = sg(u);  b = Bmax*s;
        info.b_span_ns = (max(b)-min(b))/(2*pi*p.df)*1e9;
        tau = b/(2*pi*p.df);
        tau = min(max(tau - min(tau), 0), p.tau_max);
        phi = mod(a - 2*pi*p.fc*tau, 2*pi);
        phi = reshape(phi, p.Naz, p.Nel);
        tau = reshape(tau, p.Naz, p.Nel);

    % ================= separated =====================================
    case 'sep'
        [~, ~, s0] = jp_design(p, az, el, alpha, 'sep-ls');
        Baz = Bmax/2;  Bel = Bmax/2;                 % half the range per axis
        aaz  = s0.phi_az(:) + 2*pi*p.fc*s0.tau_az(:);
        baz0 = 2*pi*p.df*s0.tau_az(:);
        ael  = s0.phi_el(:) + 2*pi*p.fc*s0.tau_el(:);
        bel0 = 2*pi*p.df*s0.tau_el(:);
        uaz = lg(min(max(baz0/Baz, 1e-3), 1-1e-3));
        uel = lg(min(max(bel0/Bel, 1e-3), 1-1e-3));

        OmA = zeros(p.Naz, Nk, Nu);  OmE = zeros(p.Nel, Nk, Nu);
        for i = 1:Nu
            aa = az(i)*pi/180;  ee = el(i)*pi/180;
            OmA(:,:,i) = (pi*(0:p.Naz-1).'*(sin(aa)*sin(ee))) * (f/p.fc);
            OmE(:,:,i) = (pi*(0:p.Nel-1).'*cos(ee))            * (f/p.fc);
        end
        MA = exp(-1j*OmA);  ME = exp(-1j*OmE);

        mAA = zeros(p.Naz,1); vAA = zeros(p.Naz,1);
        mUA = zeros(p.Naz,1); vUA = zeros(p.Naz,1);
        mEE = zeros(p.Nel,1); vEE = zeros(p.Nel,1);
        mUE = zeros(p.Nel,1); vUE = zeros(p.Nel,1);
        tt = 0;  Glhist = [];
        t0 = tic;

        for it = 1:opt.maxIter
            saz = sg(uaz);  baz = Baz*saz;
            se  = sg(uel);  bel = Bel*se;
            AA  = exp(1j*(aaz + baz*mrow)) .* MA;
            AE  = exp(1j*(ael + bel*mrow)) .* ME;
            Saz = reshape(sum(AA,1), Nk, Nu);
            Sel = reshape(sum(AE,1), Nk, Nu);
            Gk  = (abs(Saz).^2).*(abs(Sel).^2) / N;
            Gmean = sum(Gk.*selu, 1)./ns;
            Gl = sum(10*log10(Gmean));
            Glhist(end+1) = Gl;                    %#ok<AGROW>
            if converged(Glhist, opt.zeta, opt.win), break; end

            w = 10/log(10) ./ Gmean;
            QA = zeros(Nk, Nu);  QE = zeros(Nk, Nu);
            for i = 1:Nu
                k = selu(:,i);
                QA(k,i) = w(i)*(2/N)/ns(i) * conj(Saz(k,i)).*abs(Sel(k,i)).^2;
                QE(k,i) = w(i)*(2/N)/ns(i) * abs(Saz(k,i)).^2.*conj(Sel(k,i));
            end
            WA = zeros(p.Naz, Nk);  WE = zeros(p.Nel, Nk);
            for i = 1:Nu
                WA = WA + AA(:,:,i).*QA(:,i).';
                WE = WE + AE(:,:,i).*QE(:,i).';
            end
            gaaz = -imag(sum(WA, 2));   gbaz = -imag(sum(WA.*mrow, 2));
            gael = -imag(sum(WE, 2));   gbel = -imag(sum(WE.*mrow, 2));

            [aaz, mAA, vAA, tt] = adam(aaz, gaaz, mAA, vAA, tt, opt.lr, b1, b2, epsA);
            [uaz, mUA, vUA, ~ ] = adam(uaz, gbaz.*(Baz*saz.*(1-saz)), mUA, vUA, tt, opt.lr, b1, b2, epsA);
            [ael, mEE, vEE, ~ ] = adam(ael, gael, mEE, vEE, tt, opt.lr, b1, b2, epsA);
            [uel, mUE, vUE, ~ ] = adam(uel, gbel.*(Bel*se.*(1-se)), mUE, vUE, tt, opt.lr, b1, b2, epsA);

            if opt.verbose && mod(it,50) == 0
                fprintf('    gd it %3d : Gl = %.4f\n', it, Gl);
            end
        end
        info.nIter = numel(Glhist)-1;
        info.Gl_opt = Gl;

        tau_a = min(max(baz/(2*pi*p.df), 0), p.tau_max);
        tau_e = min(max(bel/(2*pi*p.df), 0), p.tau_max);
        A_tot = aaz + ael.';          % a-values only: h = A_tot + 2*pi*df*tau*m
        T_tot = tau_a + tau_e.';
        T_tot = min(max(T_tot - min(T_tot(:)), 0), p.tau_max);
        phi = mod(A_tot - 2*pi*p.fc*T_tot, 2*pi);
        tau = T_tot;
        info.b_span_ns = (max(T_tot(:))-min(T_tot(:)))*1e9;

    otherwise
        error('jp_gd: unknown mode %s', opt.mode);
end

info.time   = toc(t0);
info.Glhist = Glhist;

if opt.quantize
    qq = 2*pi/2^p.beta_ps;
    phi = mod(round(phi/qq)*qq, 2*pi);
    tau = min(max(round(tau/p.tau_p)*p.tau_p, 0), p.tau_max);
end

info.tau_span_raw = max(tau(:));
info.saturated    = max(tau(:)) > p.tau_max + eps;

R = jp_eval(p, phi, tau, az, el, blk);
info.Gl_final = R.Gl;
info.Gl_drop  = R.Gl - info.Gl_opt;      % ~0 unless the design was clipped
info.gdb      = R.gdb;
end

% =====================================================================
function tf = converged(Glhist, zeta, win)
% Convergence of Alg.1 line 18, evaluated over a window instead of a single
% step: one Adam step changes Gl by ~1e-4 relative, so the per-step test
% stops GD on its very first iteration.  Progress over the last `win`
% iterations below zeta*|Gl| is used instead (documented as gap G13).
tf = false;
if numel(Glhist) > win
    if (Glhist(end) - Glhist(end-win)) < zeta*abs(Glhist(end))
        tf = true;
    end
end
end

% ---------------------------------------------------------------------
function [x, ms, vs, t] = adam(x, g, ms, vs, t, lr, b1, b2, epsA)
t  = t + 1;
ms = b1*ms + (1-b1)*g;
vs = b2*vs + (1-b2)*g.^2;
mh = ms/(1-b1^t);
vh = vs/(1-b2^t);
x  = x + lr*mh./(sqrt(vh) + epsA);
end
