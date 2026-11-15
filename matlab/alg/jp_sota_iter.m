function [phi, tau, info] = jp_sota_iter(p, az, el, alpha, opt)
% SOTA baseline for Fig.6 / Fig.7: Algorithm 1 of
%   V. V. Ratnam, J. Mo, A. Alammouri, B. L. Ng, J. Zhang, A. F. Molisch,
%   "Joint Phase-Time Arrays: A Paradigm for Frequency-Dependent Analog
%   Beamforming in 6G", IEEE Access, vol. 10, 2022   (= this paper's ref. [5]),
% extended to the 3D array of this paper exactly as Sec. IV-B describes
% ("we have extended their iterative algorithm to 3D").
%
% ---------------------------------------------------------------------------
% WHAT [5] PRESCRIBES (eq. numbers below are [5]'s)
%
%   target      : a desired unit-norm beam bbar_k for every subcarrier k
%   objective   : max_{tau,phi,a}  sum_k omega_k * Re[ e^{j angle(a_k)} bbar_k^H T P d_k ]
%                 (eq. 4; |a_k| = ||b_k|| is fixed by eq. 5, only its phase is free)
%   alternating : for each element m (its TTD group M_n)
%                   tau_n  <- argmax_tau | sum_k omega_k e^{j angle(a_k)} [bbar_k]_m^* e^{-j2pi f_k tau} |   (8a)
%                   phi_m  <- angle( sum_k omega_k e^{j angle(a_k)} [bbar_k]_m e^{+j2pi f_k tau_n} )        (8b)
%                 for each subcarrier
%                   angle(a_k) <- angle( sum_n [bbar_k]_n e^{-jphi_m} e^{+j2pi f_k tau_n} )                (11)
%                 and (8a) may be replaced by Lemma 1: an unwrapped weighted-LS
%                 fit of the phase ramp, followed by clipping to the TTD range (9)-(10).
%
% ---------------------------------------------------------------------------
% THE 3D EXTENSION (the only degrees of freedom this paper leaves open)
%
%  (i)  one TTD per antenna element, N = Naz*Nel (this paper's hardware assumption:
%       "the number of delay elements increases linearly with the number of
%       antennas in both azimuth and elevation"), so M_n = {n} and the
%       per-element sum over m in M_n collapses to a single term.
%  (ii) the target for subcarrier k is the 3D steering vector of the user that
%       owns k:  b_k = a(f_k, theta_i(k)) / sqrt(N), i.e. exactly this paper's eq. (3)
%       array response, so that |bbar_k^H T P d_k| = sqrt(G_k) with G_k = eq. (3).
%       [5] only defines b_k for its two single-user beam behaviours; the
%       multi-user FDM target is the natural extension.
% (iii) omega_k = 1.
%
% CONSEQUENCE THAT MATTERS: with angle(a_k) free per subcarrier, the quantity
% [5] actually maximises is  sum_k sqrt(G_k)  -- an equal weight per *subcarrier*.
% this paper maximises sum_i log(mean_{k in F_i} G_k) -- an equal weight per *user*.
% At small alpha the two differ by a lot; see run_a7.m part 3.
%
% opt fields
%   variant  'lemma1' (default, eq. 9-10) | 'search' (eq. 8a line search)
%   maxIter  10    ([5] Sec. VII/VIII: "we use max-iter = 10 for Algorithm 1")
%   omega    [] (= all ones) or 1 x (M+1) weights
%   taumode  'mod'  : eq. (10a), wrap into the TTD range with mod
%            'clip' : saturate into the range instead (see note below)
%   verbose  false
%
% NOTE on (10a) vs clipping: eq.(10a) wraps tau modulo kappa/W. kappa/W is not
% a period of exp(-j2pi f_k tau) (1/df = 8.33 us is), so wrapping is a genuine
% intervention, not an equivalence. It only fires when the LS solution leaves
% the range; we run both modes and report the difference.
%
% Returns phi, tau as Naz x Nel, info with the SOTA objective history.

if nargin < 5, opt = struct(); end
opt = jp_filldef(opt, struct('variant','lemma1','maxIter',10,'omega',[], ...
                             'taumode','mod','dtgrid',[],'freePhase',true, ...
                             'verbose',false));

m  = jp_mvec(p);   f = p.fc + m*p.df;   Nk = p.M+1;
N  = p.Nelem;      Nu = numel(az);
blk = jp_block(p, alpha);
if isempty(opt.omega), om = ones(1,Nk); else, om = opt.omega(:).'; end
Tm = p.tau_max;                 % = kappa/W of [5]

% ---------------- target beamformer, one user per subcarrier --------------
Om = zeros(N, Nk);
for i = 1:Nu
    sel = (blk == i);
    nu  = jp_nu(p, az(i), el(i));  nu = nu(:);
    Om(:, sel) = nu * (f(sel)/p.fc);
end
bn = exp(1j*Om)/sqrt(N);        % N x Nk, |bn| = 1/sqrt(N)  (constant -> plain LS)
Y  = angle(bn);                 % principal-value phase, N x Nk

% ---------------- Algorithm 1 --------------------------------------------
al   = zeros(1,Nk);             % line 3 of Alg.1: initialise angle(a_k) = 0
phi  = zeros(N,1);
tau  = zeros(N,1);
F    = zeros(1,opt.maxIter+1);
nMod = 0;                       % how many elements the wrap actually moves

switch lower(opt.variant)
    case 'lemma1'
        % unknowns [phi - 2*pi*fc*tau ; b] with b = 2*pi*df*tau the phase
        % slope per subcarrier.  Scaling the second column by m (|m| <= 396)
        % instead of by 2*pi*df*m (~1e8) is essential: with the raw scaling the
        % normal matrix has condition number ~1e21 in double precision.
        X  = [ones(Nk,1), -m(:)];                % a0 - b*m  ~ unwrapped target
        Xt = X.';  XX = Xt*X;
    case 'search'
        if isempty(opt.dtgrid), dtg = p.tau_p; else, dtg = opt.dtgrid; end
        tg = (0:dtg:Tm).';                       % eq.(8a) line search grid
        Ee = exp(-1j*2*pi*(f(:)*tg(:).'));
    otherwise
        error('jp_sota_iter: unknown variant %s', opt.variant);
end

t0 = tic;
for it = 1:opt.maxIter
    % ---- eq.(9): unwrapped target phi - 2pi f_k tau_n  (N x Nk) ----------
    T = Y + al;                                  % one row per element
    T = unwrap_rows(T);

    % ---- eq.(8a) or Lemma 1 : delays ------------------------------------
    switch lower(opt.variant)
        case 'lemma1'
            coef  = XX\(Xt*T.');                 % 2 x N
            tauLS = coef(2,:).' / (2*pi*p.df);   % b = 2*pi*df*tau  ->  tau
            switch lower(opt.taumode)
                case 'mod'
                    tw = mod(tauLS, Tm);
                    nMod = nMod + nnz(abs(tw - tauLS) > 1e-12);
                    tau = tw;
                case 'clip'
                    tau = min(max(tauLS, 0), Tm);
                otherwise
                    error('jp_sota_iter: unknown taumode %s', opt.taumode);
            end
        case 'search'
            Wa = conj(bn).*exp(1j*al(:)).';      % N x Nk, eq.(8a) inner sum
            P  = Wa * Ee;                        % N x nTg
            [~,kb] = max(abs(P), [], 2);
            tau = tg(kb);
    end

    % ---- eq.(8b) : phase shifters ---------------------------------------
    % SIGN CONVENTION. [5] writes the TTD term as exp(-j2pi f_k tau_n) while
    % eq.(2) writes it as exp(+j2pi f_k tau_n).  Eqs. (7)-(11) of [5] are
    % mutually consistent only inside their own convention, so the whole
    % iteration is run in [5]'s convention (E below carries the minus) and the
    % result is mapped back at the end of the function.
    E   = exp(-1j*2*pi*(tau*f(:).'));            % N x Nk
    phi = angle( sum( exp(1j*al(:)).'.*bn.*conj(E), 2) );   % N x 1

    % ---- eq.(11) : per-subcarrier phase of the achieved beam ------------
    % opt.freePhase = false keeps angle(a_k) = 0, i.e. it drops the free
    % per-subcarrier phase of [5] eq.(5)/(11) and turns Algorithm 1 into a
    % plain alternating LS on a fixed target.  This is NOT what [5] does, but
    % it is the variant that reproduces this paper's stated 1.54 dB gap -- see run_a7.
    S    = sum( conj(bn).*E.*exp(1j*phi), 1 );   % 1 x Nk, = sqrt(G_k)
    if opt.freePhase
        al = angle(S);
    end
    F(it) = sum(om.*abs(S));

    if opt.verbose
        fprintf('    sota-%s it %2d : F = %.6f  (tau span %.1f ns)\n', ...
                opt.variant, it, F(it), (max(tau)-min(tau))*1e9);
    end
end
info.time    = toc(t0);
info.F       = F(1:opt.maxIter);
info.F_final = F(opt.maxIter);
info.nIter   = opt.maxIter;
info.variant = opt.variant;
info.freePhase = opt.freePhase;
info.taumode = opt.taumode;
info.tauLSspan_ns = (max(tau)-min(tau))*1e9;
info.nWrapped = nMod;

% describe the achieved beam with the target paper's metric (done by the caller with
% jp_eval, which recomputes G from the hardware model, so no state is hidden)
%
% Map [5]'s convention back to the target paper's: tau_paper = tau_max_5 - tau_5 turns
% exp(-j2pi f tau) into exp(+j2pi f tau') times the *common* factor
% exp(+j2pi f tau_max_5).  A factor common to all elements multiplies every
% subcarrier's array sum by a pure phase, so |G_k| is untouched, and the
% result is physically realisable (tau' >= 0).
tau = max(tau) - tau;
phi = mod(phi, 2*pi);
phi = reshape(phi, p.Naz, p.Nel);
tau = reshape(tau, p.Naz, p.Nel);
end

% =====================================================================
function T = unwrap_rows(T)
% U(.) of [5] eq.(9): add integer multiples of 2*pi along k so that adjacent
% subcarriers differ by at most pi.
d = diff(T, 1, 2);
d = d - 2*pi*round(d/(2*pi));
T = [T(:,1), T(:,1) + cumsum(d, 2)];
end
