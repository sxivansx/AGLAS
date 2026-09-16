function S = aglas_sim_setup(preset, r_command, gamma, lambda)
%AGLAS_SIM_SETUP  Build every matrix the Simulink model needs, in the base workspace.
%
%   S = AGLAS_SIM_SETUP()
%   S = AGLAS_SIM_SETUP(preset, r_command, gamma, lambda)
%
%   preset      configuration preset, default 'baseline'
%   r_command   LQR control weight, default 0.3
%   gamma       MRAC adaptation rate, default 0.5. At the design point there
%               is no model error to adapt to, so any adaptation can only add
%               control activity; 5 costs about 28 % of the achievable load
%               reduction, 0.5 costs about 4 %.
%   lambda      true control effectiveness multiplier applied to the PLANT
%               only, default 1. Use a value below 1 to simulate a degraded
%               surface the controller does not know about.
%
%   Run this before opening or simulating aglas_gla.slx. It puts a struct
%   named AGLAS in the base workspace and also assigns the individual matrices
%   the blocks reference by name, because Simulink block parameters are
%   evaluated in the base workspace.
%
%   The design matrices come from the same functions the pure-MATLAB pipeline
%   uses, so the Simulink model and CLOSED_LOOP_SIM are driven by identical
%   numbers. That is deliberate: the .m simulation is the verified reference,
%   and any disagreement between the two is a wiring error in the model rather
%   than a difference of physics.

    if nargin < 1 || isempty(preset),    preset = 'baseline'; end
    if nargin < 2 || isempty(r_command), r_command = 0.3;     end
    if nargin < 3 || isempty(gamma),     gamma = 0.5;         end
    if nargin < 4 || isempty(lambda),    lambda = 1;          end

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(here), 'src', 'common'));

    cfg = aglas_config(preset);
    cfg.lqr.r_command = r_command;

    fem   = beam_fem(cfg);
    modes = modal_analysis(fem, cfg.fem.n_modes);
    cs    = control_surface(cfg, fem, modes);
    P     = plant_statespace(cfg, fem, modes, cs);

    base = lqg_design(cfg, P);
    mrac = mrac_design(cfg, P, base, struct('gamma', gamma, 'sigma', 0.05));

    ns = P.n_states;
    na = base.n_aug;

    % --- true plant, with the unknown effectiveness applied ---------------
    Bt = P.B * lambda;

    % Single State-Space block: inputs [u; w], outputs [y(2); z(2)]
    S.plant.A  = P.A;
    S.plant.B  = [Bt, P.Bw];
    S.plant.C  = [P.C; P.Cz];
    S.plant.D  = [P.D, P.Dw; P.Dz, P.Dzw];
    S.plant.x0 = zeros(ns, 1);

    % --- estimator, written as a plain state-space block -------------------
    % xhat_dot = (Ae - L*Ce)*xhat + (Be - L*De)*u + L*y
    S.est.A  = base.Ae - base.L*base.Ce;
    S.est.B  = [base.Be - base.L*base.De, base.L];
    S.est.C  = eye(na);
    S.est.D  = zeros(na, 3);
    S.est.x0 = zeros(na, 1);

    % --- reference model for MRAC ------------------------------------------
    S.ref.A  = mrac.Am;
    S.ref.B  = P.Bw;
    S.ref.C  = eye(ns);
    S.ref.D  = zeros(ns, 1);
    S.ref.x0 = zeros(ns, 1);

    % --- gains and limits ---------------------------------------------------
    S.K          = base.K;
    S.n_states   = ns;
    S.n_aug      = na;
    S.idx_gust   = base.idx_gust;
    S.delta_max  = deg2rad(cfg.control.delta_max);
    S.rate_max   = deg2rad(cfg.control.rate_max);
    S.Gamma      = mrac.Gamma;
    S.sigma      = mrac.sigma;
    S.theta_max  = mrac.theta_max;
    S.PB         = mrac.PB;
    S.normalise  = mrac.normalise;
    S.lambda     = lambda;

    % --- gust signal --------------------------------------------------------
    dt = 2e-4;
    t  = (0:dt:6).';
    wg = gust_one_minus_cos(t.', cfg.gust.U_ds, cfg.gust.t_g).';
    S.gust    = [t, wg];
    S.dt      = dt;
    S.t_final = t(end);

    % 1 enables the adaptive branch, 0 leaves plain LQG. Kept as a gain rather
    % than a rebuild so the two can be compared without regenerating the model.
    S.adapt_on = 1;

    S.cfg = cfg; S.P = P; S.base = base; S.mrac = mrac;

    % Simulink evaluates block parameters in the base workspace.
    assignin('base', 'AGLAS', S);
    fprintf('AGLAS Simulink parameters ready in the base workspace.\n');
    fprintf('  plant  %d states, effectiveness lambda = %.2f\n', ns, lambda);
    fprintf('  estimator %d states (structure + gust)\n', na);
    fprintf('  actuator limits %.1f deg, %.0f deg/s\n', ...
            cfg.control.delta_max, cfg.control.rate_max);
    fprintf('  simulate to t = %.1f s with a fixed step of %.0e s\n', S.t_final, dt);
end
