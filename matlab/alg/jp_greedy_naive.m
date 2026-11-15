function [phi, tau, info] = jp_greedy_naive(p, az, el, alpha, opt)
% Literal implementation of Algorithm 1, used ONLY to verify that the
% incremental-evaluation version (jp_greedy.m) is numerically identical.
% Every candidate is scored by a full recomputation of eq. (5) through
% jp_eval, so the cost per candidate is O(N_elem) exactly as in a direct
% reading of the algorithm.  Do not use it for timing comparisons with GD
% unless that is the point being made.
%
%   opt.mode  'joint' | 'sep' ;  opt.zeta ; opt.maxIter

if nargin < 5, opt = struct(); end
opt = jp_filldef(opt, struct('mode','joint','zeta',1e-4,'maxIter',10,'verbose',false));

blk = jp_block(p, alpha);
q   = 2*pi/2^p.beta_ps;
tg  = (0:p.tau_p:p.tau_max).';
pg  = (0:q:(2*pi-q)).';

switch lower(opt.mode)
    case 'joint'
        [phi, tau] = jp_design(p, az, el, alpha, 'joint-ls');
        tau = min(max(round(tau/p.tau_p)*p.tau_p, 0), p.tau_max);
        phi = mod(round(phi/q)*q, 2*pi);

        Gl = evalGl(phi, tau);
        t0 = tic;
        for it = 1:opt.maxIter
            Gl0 = Gl;
            for e = 1:p.Nelem
                best = -Inf;  barg = tau(e);
                for c = 1:numel(tg)
                    tau(e) = tg(c);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = tg(c); end
                end
                tau(e) = barg;
                best = -Inf;  barg = phi(e);
                for c = 1:numel(pg)
                    phi(e) = pg(c);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = pg(c); end
                end
                phi(e) = barg;
            end
            Gl = evalGl(phi, tau);
            if abs(Gl - Gl0) < opt.zeta*abs(Gl), break; end
            if opt.verbose
                fprintf('    naive it %2d : Gl = %.4f (d = %+.4f)\n', it, Gl, Gl-Gl0);
            end
        end

    case 'sep'
        [~, ~, s] = jp_design(p, az, el, alpha, 'sep-ls');
        phi_a = mod(round(s.phi_az(:)/q)*q, 2*pi);
        tau_a = min(max(round(s.tau_az(:)/p.tau_p)*p.tau_p, 0), p.tau_max);
        phi_e = mod(round(s.phi_el(:)/q)*q, 2*pi);
        tau_e = min(max(round(s.tau_el(:)/p.tau_p)*p.tau_p, 0), p.tau_max);
        asm = @(pa,ta,pe,te) deal(mod(pa+pe.', 2*pi), min(ta+te.', p.tau_max));
        [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);

        Gl = evalGl(phi, tau);
        t0 = tic;
        for it = 1:opt.maxIter
            Gl0 = Gl;
            for e = 1:p.Naz
                best = -Inf;  barg = tau_a(e);
                for c = 1:numel(tg)
                    tau_a(e) = tg(c);
                    if max(tau_a)+max(tau_e) > p.tau_max + 1e-15, continue; end
                    [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = tg(c); end
                end
                tau_a(e) = barg;
                best = -Inf;  barg = phi_a(e);
                for c = 1:numel(pg)
                    phi_a(e) = pg(c);
                    [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = pg(c); end
                end
                phi_a(e) = barg;
            end
            for e = 1:p.Nel
                best = -Inf;  barg = tau_e(e);
                for c = 1:numel(tg)
                    tau_e(e) = tg(c);
                    if max(tau_a)+max(tau_e) > p.tau_max + 1e-15, continue; end
                    [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = tg(c); end
                end
                tau_e(e) = barg;
                best = -Inf;  barg = phi_e(e);
                for c = 1:numel(pg)
                    phi_e(e) = pg(c);
                    [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);
                    g = evalGl(phi, tau);
                    if g > best, best = g;  barg = pg(c); end
                end
                phi_e(e) = barg;
            end
            [phi, tau] = asm(phi_a, tau_a, phi_e, tau_e);
            Gl = evalGl(phi, tau);
            if abs(Gl - Gl0) < opt.zeta*abs(Gl), break; end
            if opt.verbose
                fprintf('    naive it %2d : Gl = %.4f (d = %+.4f)\n', it, Gl, Gl-Gl0);
            end
        end

    otherwise
        error('jp_greedy_naive: unknown mode %s', opt.mode);
end

info.mode   = opt.mode;
info.Gl     = Gl;
info.nIter  = it;
info.time   = toc(t0);
info.tau_span_raw = max(tau(:));

    function g = evalGl(ph, ta)
        R = jp_eval(p, ph, ta, az, el, blk);
        g = R.Gl;
    end
end
