function compare_simulink(mdl)
%COMPARE_SIMULINK  Check the Simulink model against the verified .m reference.
%
%   compare_simulink
%   COMPARE_SIMULINK(name)
%
%   Runs the Simulink model and CLOSED_LOOP_SIM on identical inputs and reports
%   the difference. Both are driven from the same AGLAS struct, so the physics
%   is identical by construction and any disagreement is a wiring or
%   configuration error in the model.
%
%   This exists because a block diagram is easy to get subtly wrong: a swapped
%   Mux order, a Selector off by one, a sign on a Sum. None of those raise an
%   error, they just give a plausible wrong answer. The .m path is covered by
%   the test suite, so it is the reference and the model is the thing under
%   test, not the other way round.
%
%   One difference between the two is structural and will not go away. The .m
%   simulation computes the control once per step and holds it across the four
%   RK4 stages, which is what a digital controller sampling at 5 kHz actually
%   does. Simulink integrates the controller as part of one continuous system,
%   so its command varies within the step. Expect the two to agree closely on
%   peaks and RMS and to drift by a few percent instantaneously wherever the
%   response is lightly damped and the signals are moving fastest. That is a
%   modelling difference, not an error in either.
%
%   Diagnostic tip: set AGLAS.adapt_on = 0 and re-run. That removes the two
%   adaptive integrators from the loop entirely. If plain LQG agrees tightly
%   and MRAC does not, the discrepancy is in the adaptive path rather than in
%   the plant or estimator wiring.

    if nargin < 1 || isempty(mdl), mdl = 'aglas_gla'; end

    S     = evalin('base', 'AGLAS');
    paths = aglas_paths();

    fprintf('Running Simulink model %s ...\n', mdl);
    simOut = sim(mdl);

    m_sl = get_logged(simOut, 'aglas_moment');
    d_sl = get_logged(simOut, 'aglas_delta_cmd');

    t_sl   = m_sl.time(:).';
    mom_sl = squeeze(m_sl.signals.values).';
    del_sl = squeeze(d_sl.signals.values).';

    fprintf('Running the .m reference ...\n');
    P = S.P;
    P.B = P.B * S.lambda;
    if S.adapt_on
        ctrl = S.mrac;
    else
        ctrl = S.base;
    end
    t_m  = S.gust(:,1).';
    wg_m = S.gust(:,2).';
    ref  = closed_loop_sim(S.cfg, P, ctrl, t_m, wg_m);

    mom_ref = interp1(ref.t, ref.moment,    t_sl, 'linear', 'extrap');
    del_ref = interp1(ref.t, ref.delta_cmd, t_sl, 'linear', 'extrap');

    % Three separate measures. A single max-pointwise number is misleading:
    % two traces can agree on the peak load, which is the engineering answer,
    % while differing during a lightly damped secondary oscillation where a
    % small timing shift produces a large instantaneous difference.
    pk_mom = abs(max(abs(mom_sl)) - max(abs(mom_ref))) / max(abs(mom_ref));
    pk_del = abs(max(abs(del_sl)) - max(abs(del_ref))) / max(max(abs(del_ref)), eps);
    rms_mom = sqrt(mean((mom_sl - mom_ref).^2)) / sqrt(mean(mom_ref.^2));
    mx_mom  = max(abs(mom_sl - mom_ref)) / max(abs(mom_ref));

    fprintf('\n%-30s %13s %13s\n', 'quantity', 'Simulink', 'reference');
    fprintf('%-30s %13.2f %13.2f\n', 'peak root moment [kN m]', ...
            max(abs(mom_sl))/1e3, max(abs(mom_ref))/1e3);
    fprintf('%-30s %13.2f %13.2f\n', 'peak command [deg]', ...
            rad2deg(max(abs(del_sl))), rad2deg(max(abs(del_ref))));

    fprintf('\n%-30s %9.3f %%\n', 'peak moment difference',   100*pk_mom);
    fprintf('%-30s %9.3f %%\n',   'peak command difference',  100*pk_del);
    fprintf('%-30s %9.3f %%\n',   'waveform RMS difference',  100*rms_mom);
    fprintf('%-30s %9.3f %%\n',   'worst pointwise difference', 100*mx_mom);

    tol_peak = 0.02;
    tol_rms  = 0.05;
    ok = (pk_mom < tol_peak) && (pk_del < tol_peak) && (rms_mom < tol_rms);

    fprintf('\n');
    if ok
        fprintf('Model agrees with the reference. Peaks within %.0f %%, waveform within %.0f %% RMS.\n', ...
                100*tol_peak, 100*tol_rms);
        if mx_mom > 0.05
            fprintf(['Worst pointwise difference is %.1f %%, concentrated where the\n' ...
                     'response is lightly damped. A small timing shift there is\n' ...
                     'expected and does not affect the load the design is sized by.\n'], ...
                     100*mx_mom);
        end
    else
        fprintf('DISAGREEMENT. Check, in this order:\n');
        fprintf('  1. AdaptIn Mux order, must be [xhat(9); xm(8); theta(8)]\n');
        fprintf('  2. UadIn Mux order, [theta(8); xhat_struct(8)]\n');
        fprintf('  3. PickStruct indices 1:8, PickGust index 9\n');
        fprintf('  4. SumU signs, both +\n');
        fprintf('  5. Kgain Multiplication set to Matrix(K*u)\n');
        fprintf('  6. Saturation placed BEFORE Rate Limiter, not after\n');
    end

    % ---------------------------------------------------------------- figure
    fig = figure('Name', 'Simulink against reference', 'Color', 'w', ...
                 'Position', [100 100 900 620]);
    subplot(2,1,1);
    plot(t_sl, mom_sl/1e3, 'LineWidth', 1.8); hold on;
    plot(t_sl, mom_ref/1e3, '--', 'LineWidth', 1.5); grid on;
    ylabel('root moment [kN m]');
    legend('Simulink', '.m reference', 'Location', 'northeast');
    title(sprintf('Model verification: %s, adaptation %s, lambda = %.2f', ...
          mdl, ternary(S.adapt_on ~= 0, 'on', 'off'), S.lambda), ...
          'Interpreter', 'none');

    subplot(2,1,2);
    plot(t_sl, rad2deg(del_sl), 'LineWidth', 1.8); hold on;
    plot(t_sl, rad2deg(del_ref), '--', 'LineWidth', 1.5);
    plot(t_sl([1 end]),  S.cfg.control.delta_max*[1 1], 'k:', 'LineWidth', 1.2);
    plot(t_sl([1 end]), -S.cfg.control.delta_max*[1 1], 'k:', 'LineWidth', 1.2);
    grid on; xlabel('time [s]'); ylabel('command [deg]');
    legend('Simulink', '.m reference', 'deflection limit', 'Location', 'southeast');

    fig_file = fullfile(paths.results, 'simulink_verification.png');
    print(fig, fig_file, '-dpng', '-r150');

    % ------------------------------------------------------------------ data
    sl.t          = t_sl;
    sl.moment     = mom_sl;
    sl.delta_cmd  = del_sl;
    sl.tip        = try_logged(simOut, 'aglas_tip');
    sl.theta      = try_logged(simOut, 'aglas_theta');
    sl.u_adaptive = try_logged(simOut, 'aglas_u_ad');
    sl.gust       = try_logged(simOut, 'aglas_gust');

    comparison.peak_moment_pct  = 100*pk_mom;
    comparison.peak_command_pct = 100*pk_del;
    comparison.rms_pct          = 100*rms_mom;
    comparison.max_pointwise_pct= 100*mx_mom;
    comparison.agrees           = ok;

    settings.adapt_on = S.adapt_on;
    settings.lambda   = S.lambda;
    settings.r_command= S.cfg.lqr.r_command;
    settings.dt       = S.dt;

    data_file = fullfile(paths.data, 'simulink_run.mat');
    reference = ref; %#ok<NASGU>
    save(data_file, 'sl', 'reference', 'comparison', 'settings');

    fprintf('\nSaved %s\n', fig_file);
    fprintf('Saved %s\n', data_file);
end

% ------------------------------------------------------------------------
function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

% ------------------------------------------------------------------------
function v = try_logged(simOut, name)
%TRY_LOGGED  Fetch an optional signal, returning empty rather than erroring.
    v = [];
    try %#ok<TRYNC>
        v = get_logged(simOut, name);
    end
end

% ------------------------------------------------------------------------
function sig = get_logged(simOut, name)
%GET_LOGGED  Retrieve a To Workspace signal, whichever way this release returns it.
%
%   Recent MATLAB returns a Simulink.SimulationOutput object; older
%   configurations drop the variables straight into the base workspace. Try the
%   object first, fall back to the base workspace, and fail with a message that
%   names the missing signal rather than letting the lookup resolve to some
%   unrelated built-in function of the same name.

    sig = [];

    try %#ok<TRYNC>
        if isa(simOut, 'Simulink.SimulationOutput')
            avail = simOut.who;
            if any(strcmp(avail, name))
                sig = simOut.get(name);
            end
        end
    end

    if isempty(sig)
        try %#ok<TRYNC>
            if evalin('base', sprintf('exist(''%s'', ''var'')', name)) == 1
                sig = evalin('base', name);
            end
        end
    end

    if isempty(sig) || ~isstruct(sig) || ~isfield(sig, 'time')
        error('compare_simulink:missingLog', ...
            ['The model did not log "%s". Check that the To Workspace block ' ...
             'exists, that its Variable name is "%s", and that its Save ' ...
             'format is "Structure With Time". Regenerating with ' ...
             'build_gla_model will restore all six log blocks.'], name, name);
    end
end
