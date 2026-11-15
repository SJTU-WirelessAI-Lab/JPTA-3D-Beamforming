function [phi, tau, info] = jp_design(p, az, el, alpha, mode)
% Solve the JPTA beam design problem of eq. (6) analytically.
%
% mode:
%   'joint-ls'  : joint phase/delay, least squares      (Sec. III-A-1, eq. 11-13)
%   'joint-inf' : joint phase/delay, infinity norm      (Sec. III-A-2)
%   'sep-ls'    : separated azimuth/elevation, LS       (Sec. III-B)
%   'sep-inf'   : separated azimuth/elevation, inf norm
%
% Returns quantized PS phases phi (Naz x Nel, rad in [0,2pi)) and TTD delays
% tau (Naz x Nel, s), plus diagnostics in info.

p.M;                                    %#ok<*NASGU>  (params passed whole)
m   = jp_mvec(p);
blk = jp_block(p, alpha);
Nu  = numel(az);

switch mode
    case {'joint-ls','joint-inf'}
        % ---- Sec. III-A : one equation per antenna element --------
        if isfield(p,'nowrap') && p.nowrap
            K = zeros(p.Naz, p.Nel, Nu);            % no 2*pi*K unwrapping
        else
            K = jp_unwrap_k(p, az, el);             % Naz x Nel x Nu
        end
        nu_list = cell(1,Nu);
        for i = 1:Nu
            nu_list{i} = jp_nu(p, az(i), el(i)) + 2*pi*K(:,:,i);
        end
        B = jp_target_blocks(p, nu_list, blk);      % (M+1) x Nelem
        if strcmp(mode,'joint-ls')
            [x1, x2] = jp_ls_fit(m, B);
        else
            [x1, x2] = jp_inf_fit(m, B);
        end
        [phi, tau, info] = jp_finalize(p, x1, x2, [p.Naz p.Nel]);
        info.x1 = x1;  info.x2 = x2;                % raw LS solution

    case {'sep-ls','sep-inf'}
        % ---- Sec. III-B : azimuth and elevation designed separately
        nu_az = cell(1,Nu);  nu_el = cell(1,Nu);
        for i = 1:Nu
            a = az(i)*pi/180;  e = el(i)*pi/180;
            nu_az{i} = pi*(0:p.Naz-1).' * (sin(a)*sin(e));   % Naz x 1
            nu_el{i} = pi*(0:p.Nel-1).' * cos(e);            % Nel x 1
        end
        if isfield(p,'nowrap') && p.nowrap
            Kaz = zeros(p.Naz, Nu);  Kel = zeros(p.Nel, Nu);
        else
            Kaz = jp_unwrap_k_1d(nu_az);
            Kel = jp_unwrap_k_1d(nu_el);
        end
        for i = 1:Nu
            nu_az{i} = nu_az{i} + 2*pi*Kaz(:,i);
            nu_el{i} = nu_el{i} + 2*pi*Kel(:,i);
        end
        Baz = jp_target_blocks(p, nu_az, blk);      % (M+1) x Naz
        Bel = jp_target_blocks(p, nu_el, blk);      % (M+1) x Nel

        if strcmp(mode,'sep-ls')
            [x1a, x2a] = jp_ls_fit(m, Baz);
            [x1e, x2e] = jp_ls_fit(m, Bel);
        else
            [x1a, x2a] = jp_inf_fit(m, Baz);
            [x1e, x2e] = jp_inf_fit(m, Bel);
        end
        [phia, taua, ia] = jp_finalize(p, x1a, x2a, [p.Naz 1]);
        [phie, taue, ie] = jp_finalize(p, x1e, x2e, [p.Nel 1]);

        phi = phia + phie.';                        % Naz x Nel
        tau = taua + taue.';                        % Naz x Nel
        phi = mod(phi, 2*pi);
        info.tau_span_raw = max(tau(:));
        info.saturated    = max(tau(:)) > p.tau_max + eps;
        if info.saturated
            tau = min(tau, p.tau_max);
        end
        info.span_az = ia.tau_span_raw;
        info.span_el = ie.tau_span_raw;
        % per-axis quantized components, needed to initialise greedy / GD
        info.phi_az = phia;  info.tau_az = taua;
        info.phi_el = phie;  info.tau_el = taue;
        info.x_az = [x1a; x2a];  info.x_el = [x1e; x2e];

    otherwise
        error('unknown mode %s', mode);
end
end
