% Copyright © 2023 Martin Schonger
% Copyright © 2025 Simone Silenzi
% This software is licensed under the GPLv3.
function options = fvbsettings(varargin)
% Collect options for the f/V/B problem. Minimal, compatible with SISBMI path.

p = inputParser; p.KeepUnmatched = true;

addParameter(p, 'seed', 13);

addParameter(p, 'enable_barrier', true);
addParameter(p, 'epsilon', 1e-4);
addParameter(p, 'init_Bc_var', true);
addParameter(p, 'constraint_version', 1);   % 1: hard; 2: slack*(eps-B); 3: slack*(eps-B^2)
addParameter(p, 'enable_regularization', true);

addParameter(p, 'deg_f', 2);
addParameter(p, 'deg_V', 2);
addParameter(p, 'deg_B', 2);
addParameter(p, 'deg_B_slack', 2);

addParameter(p, 'enable_extra_constraint', true);
addParameter(p, 'regularization_factor', 0.01);

addParameter(p, 'dataset', 'lasa');
addParameter(p, 'dataset_opts', struct('idx', 14));

% New: YALMIP options (default MOSEK). Legacy kept for compatibility (unused).
addParameter(p, 'sdpoptions', sdpsettings('solver','mosek','verbose',1,'cachesolvers',1));
addParameter(p, 'sdpoptions_penbmi', struct());

parse(p, varargin{:});
options = p.Results;

% Propagate unmatched things (e.g., SIS tuning, restrict_to_convex, keep_fc)
if ~isempty(fieldnames(p.Unmatched))
    options.unmatched = p.Unmatched;
end
end
