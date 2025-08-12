% Copyright © 2023 Martin Schonger
% Copyright © 2025 Simone Silenzi
% This software is licensed under the GPLv3.
function [f_fh, V_fh, dVdx_fh, B_fh, dBdx_fh, debug_output, fc, Vc, Bc] = fvb(rd, restrict_to_convex, xi, initial_set, unsafe_set, options)
% Solve for a polynomial dynamical system f, Lyapunov function V, and barrier B.
% SISBMI (YALMIP+MOSEK) replaces BMI/PENBMI when restrict_to_convex == 1.
%
% Inputs:
%   rd: reference data with fields M (#states), T (#samples), Data([x; xdot])
%   restrict_to_convex: 0 => convex init with fixed V = xi'*xi, 1 => SISBMI alternating
%   xi: sdpvar column vector of states (size Mx1)
%   initial_set, unsafe_set: cell arrays of polynomials >= 0 (can be empty)
%   options: struct from fvbsettings (expects options.sdpoptions for YALMIP)
%
% Outputs: function handles for f, V, dVdx, B, dBdx and coefficient vectors.

epsilon = options.epsilon;
debug_output = struct;

% Dimensions
M = rd.M;
T = rd.T;

%% Decision variables for f
deg_f = options.deg_f;
f = []; fc_var = [];
for m = 1:M
    [f_m, fc_m] = polynomial(xi, deg_f, 1);
    f      = [f; f_m];
    fc_var = [fc_var; fc_m];
end

% Monomial basis (for warm-start plumbing used elsewhere)
try
    debug_output.f_monomials = monolist(xi, deg_f);
catch
    debug_output.f_monomials = [];
end

%% Lyapunov V
if restrict_to_convex == 0
    V = xi' * xi;               % fixed quadratic for convex warm-start
    [Vc_var, ~] = coefficients(V, xi); % keep variable in scope
else
    [V, Vc_var] = polynomial(xi, options.deg_V, 1);
end
dVdx = jacobian(V, xi)';

%% Barrier B
[B, Bc_var] = polynomial(xi, options.deg_B);
dBdx = jacobian(B, xi)';

%% Common fitting pieces
xi_dot = sdpvar(M, T, 'full');
for t = 1:T
    xi_dot(:, t) = replace(f, xi, rd.Data(1:M, t));
end
xi_dot_error = xi_dot - rd.Data(M+1:end, :);
mse = sum(sum(xi_dot_error.^2)) / (2 * T);
Objective_fit = mse;

deg_B_slack = options.deg_B_slack; % still used for some slack polynomials
epxi = epsilon * sum(xi.^2, 1);

% Normalize empty sets
if nargin < 4 || isempty(initial_set), initial_set = {}; end
if nargin < 5 || isempty(unsafe_set),  unsafe_set  = {}; end

%% Domain polynomial gX >= 0 (for S-procedure multipliers)
% Use provided domain if present; otherwise build a ball covering the data.
domain_set = {};
if isfield(options,'unmatched') && isfield(options.unmatched,'domain_set') && ~isempty(options.unmatched.domain_set)
    % Expect a cell array of polynomials g_i(x) >= 0
    domain_set = options.unmatched.domain_set;
else
    X = rd.Data(1:M, :);
    R2 = 1.1 * max(sum(X.^2, 1));   % margin
    gX = R2 - sum(xi.^2);           % ball: ||x||^2 <= R^2  <=> gX >= 0
    domain_set = {gX};
end
debug_output.domain_set = domain_set;

% Helper: degree-balanced domain multipliers for SOS(poly - sum mu_i*g_i)
    function [ConOut, mults] = add_domain_sos_auto(ConIn, poly)
        ConOut = ConIn;
        mults  = {};
        if isempty(domain_set)
            ConOut = [ConOut, sos(poly)];
            return
        end
        d_poly = degree(poly, xi);
        s = 0;
        for ii = 1:numel(domain_set)
            d_g = degree(domain_set{ii}, xi);
            d_mu = max(0, d_poly - d_g);             % balance degrees so mu*g can cancel leading terms
            [mu_ii, ~] = polynomial(xi, d_mu);
            ConOut = [ConOut, sos(mu_ii)];
            s = s + mu_ii * domain_set{ii};
            mults{end+1} = mu_ii;
        end
        ConOut = [ConOut, sos(poly - s)];
    end

%% SISBMI (full alternating SOS) or convex init
if restrict_to_convex == 1
    % SIS tuning
    sis_maxit = 20;
    sis_tol   = 1e-3;
    if isfield(options, 'unmatched')
        if isfield(options.unmatched, 'sis_max_iterations'), sis_maxit = options.unmatched.sis_max_iterations; end
        if isfield(options.unmatched, 'sis_tolerance'),     sis_tol   = options.unmatched.sis_tolerance;     end
    end
    if isfield(options, 'sdpoptions') && ~isempty(options.sdpoptions)
        sdp_options = options.sdpoptions;
    else
        sdp_options = sdpsettings('solver','mosek','verbose',1,'cachesolvers',1);
    end

    % Warm start f
    fc_k = value(fc_var);
    if isempty(fc_k) || all(fc_k == 0), fc_k = zeros(size(fc_var)); end
    best_mse = inf;

    debug_output.sis = struct('it',[],'mse',[],'infoV',{{}},'infof',{{}});

    for it = 1:sis_maxit
        % ===== (A) VB-step: fix f := f_k, solve V,B =====
        f_fixed = replace(f, fc_var, fc_k);
        Constraints_VB = [];

        % Lyapunov positivity and decay (decay relaxed to domain with degree balance)
        Constraints_VB = [Constraints_VB, sos(V - epxi)];
        Vdot_fixed = sum(jacobian(V, xi)'.*f_fixed, 1);
        [Constraints_VB, ~] = add_domain_sos_auto(Constraints_VB, -Vdot_fixed - epxi);

        if options.enable_barrier
            % Initial set: B <= 0 on {g >= 0}  ->  -B - sum tau*g is SOS
            if ~isempty(initial_set)
                sos_safe = -B;
                for p = 1:length(initial_set)
                    [tau_p, ~] = polynomial(xi, deg_B_slack);
                    Constraints_VB = [Constraints_VB, sos(tau_p)];
                    sos_safe = sos_safe - tau_p * initial_set{p};
                end
                Constraints_VB = [Constraints_VB, sos(sos_safe)];
            end

            % Unsafe set: B >= epsilon on {h >= 0} ->  B - epsilon - sum sigma*h is SOS
            if ~isempty(unsafe_set)
                sos_unsafe = B - epsilon;
                for q = 1:length(unsafe_set)
                    [sig_q, ~] = polynomial(xi, deg_B_slack);
                    Constraints_VB = [Constraints_VB, sos(sig_q)];
                    sos_unsafe = sos_unsafe - sig_q * unsafe_set{q};
                end
                Constraints_VB = [Constraints_VB, sos(sos_unsafe)];
            end

            % Barrier derivative (relaxed to domain with degree balance)
            Bdot_fixed = sum(jacobian(B, xi)'.*f_fixed, 1);
            switch options.constraint_version
                case 1
                    [Constraints_VB, ~] = add_domain_sos_auto(Constraints_VB, -Bdot_fixed - epxi);
                case 2
                    [slackvar, ~] = polynomial(xi, deg_B_slack);
                    Constraints_VB = [Constraints_VB, sos(slackvar)];
                    [Constraints_VB, ~] = add_domain_sos_auto(Constraints_VB, -Bdot_fixed - epxi - slackvar*(epsilon - B));
                case 3
                    [slackvar, ~] = polynomial(xi, deg_B_slack);
                    Constraints_VB = [Constraints_VB, sos(slackvar)];
                    [Constraints_VB, ~] = add_domain_sos_auto(Constraints_VB, -Bdot_fixed - epxi - slackvar*(epsilon - B.^2));
                otherwise
                    [Constraints_VB, ~] = add_domain_sos_auto(Constraints_VB, -Bdot_fixed - epxi);
            end
        end

        [sol_VB, ~, Q_VB] = solvesos(Constraints_VB, 0, sdp_options, [Vc_var(:); Bc_var(:)]);
        debug_output.sis.it(end+1)    = it;
        debug_output.sis.infoV{end+1} = sol_VB;

        % Freeze V,B
        if isa(Vc_var,'sdpvar') && ~isempty(Vc_var)
            Vk = replace(V, Vc_var, value(Vc_var));
        else
            Vk = V; % fixed V
        end
        Bk    = replace(B, Bc_var, value(Bc_var));
        dVdxk = jacobian(Vk, xi)';
        dBdxk = jacobian(Bk, xi)';

        % ===== (B) f-step: fix V:=Vk, B:=Bk, solve f =====
        Constraints_f = [];

        Vdot_k = sum(dVdxk.*f, 1);
        [Constraints_f, ~] = add_domain_sos_auto(Constraints_f, -Vdot_k - epxi);

        if options.enable_barrier
            Bdot_k = sum(dBdxk.*f, 1);
            switch options.constraint_version
                case 1
                    [Constraints_f, ~] = add_domain_sos_auto(Constraints_f, -Bdot_k - epxi);
                case 2
                    [slackvar_f, ~] = polynomial(xi, deg_B_slack);
                    Constraints_f = [Constraints_f, sos(slackvar_f)];
                    [Constraints_f, ~] = add_domain_sos_auto(Constraints_f, -Bdot_k - epxi - slackvar_f*(epsilon - Bk));
                case 3
                    [slackvar_f, ~] = polynomial(xi, deg_B_slack);
                    Constraints_f = [Constraints_f, sos(slackvar_f)];
                    [Constraints_f, ~] = add_domain_sos_auto(Constraints_f, -Bdot_k - epxi - slackvar_f*(epsilon - Bk.^2));
                otherwise
                    [Constraints_f, ~] = add_domain_sos_auto(Constraints_f, -Bdot_k - epxi);
            end
        end

        Objective_f = Objective_fit;
        if options.enable_regularization
            Objective_f = Objective_f + options.regularization_factor * norm(fc_var(:), 2)^2;
        end

        [sol_f, ~, Q_f] = solvesos(Constraints_f, Objective_f, sdp_options, fc_var(:));
        mse_k = value(mse);
        debug_output.sis.mse(end+1) = mse_k;

        if abs(best_mse - mse_k) < sis_tol
            break;
        end
        best_mse = mse_k;
        fc_k = value(fc_var);
    end

    debug_output.Q = {Q_VB};
else
    % ===== Convex init (fixed V = xi'xi) =====
    if isfield(options,'sdpoptions') && ~isempty(options.sdpoptions)
        sdp_options = options.sdpoptions;
    else
        sdp_options = sdpsettings('solver','mosek','verbose',1,'cachesolvers',1);
    end

    Constraints = [];
    Constraints = [Constraints, sos(V - epxi)];
    Vdot = sum(dVdx.*f, 1);
    [Constraints, ~] = add_domain_sos_auto(Constraints, -Vdot - epxi);

    if options.enable_barrier
        if ~isempty(initial_set)
            sos_safe = -B;
            for p = 1:length(initial_set)
                [tau_p, ~] = polynomial(xi, deg_B_slack);
                Constraints = [Constraints, sos(tau_p)];
                sos_safe = sos_safe - tau_p * initial_set{p};
            end
            Constraints = [Constraints, sos(sos_safe)];
        end

        if ~isempty(unsafe_set)
            sos_unsafe = B - epsilon;
            for q = 1:length(unsafe_set)
                [sig_q, ~] = polynomial(xi, deg_B_slack);
                Constraints = [Constraints, sos(sig_q)];
                sos_unsafe = sos_unsafe - sig_q * unsafe_set{q};
            end
            Constraints = [Constraints, sos(sos_unsafe)];
        end

        Bdot = sum(dBdx.*f, 1);
        switch options.constraint_version
            case 1
                [Constraints, ~] = add_domain_sos_auto(Constraints, -Bdot - epxi);
            case 2
                [slackvar, ~] = polynomial(xi, deg_B_slack);
                Constraints = [Constraints, sos(slackvar)];
                [Constraints, ~] = add_domain_sos_auto(Constraints, -Bdot - epxi - slackvar*(epsilon - B));
            case 3
                [slackvar, ~] = polynomial(xi, deg_B_slack);
                Constraints = [Constraints, sos(slackvar)];
                [Constraints, ~] = add_domain_sos_auto(Constraints, -Bdot - epxi - slackvar*(epsilon - B.^2));
            otherwise
                [Constraints, ~] = add_domain_sos_auto(Constraints, -Bdot - epxi);
        end
    end

    Objective = Objective_fit;
    if options.enable_regularization
        Objective = Objective + options.regularization_factor * norm(fc_var(:), 2)^2;
    end

    [sol, ~, Q] = solvesos(Constraints, Objective, sdp_options, fc_var(:));
    debug_output.convex_init = sol;
end

%% Outputs: function handles and coefficient values
fc = value(fc_var);
f_num = replace(f, fc_var, fc);

% V may be fixed (no decision variables) in convex-init
Vc = [];
if exist('Vc_var','var') && isa(Vc_var,'sdpvar') && ~isempty(Vc_var)
    Vc    = value(Vc_var);
    V_num = replace(V, Vc_var, Vc);
else
    V_num = V;  % fixed V = xi'*xi
end

% Barrier is always parametric here
Bc    = value(Bc_var);
B_num = replace(B, Bc_var, Bc);

dVdx_num = jacobian(V_num, xi)';
dBdx_num = jacobian(B_num, xi)';

% Function handles (relies on sdpvar2fun existing in repo)
f_fhm = cell(M,1);
for m = 1:M
    f_fhm{m} = sdpvar2fun(f_num(m), xi);
end
f_fh    = @(x) cell2mat(cellfun(@(g) g(x), f_fhm, 'UniformOutput', false));
V_fh    = sdpvar2fun(V_num, xi);
dVdx_fh = @(x) double(replace(dVdx_num, xi, x));
B_fh    = sdpvar2fun(B_num, xi);
dBdx_fh = @(x) double(replace(dBdx_num, xi, x));
end
