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

    if nargin < 1 || isempty(mdl), mdl = 'aglas_gla'; end

    S = evalin('base', 'AGLAS');

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

    e_mom = max(abs(mom_sl - mom_ref)) / max(abs(mom_ref));
    e_del = max(abs(del_sl - del_ref)) / max(max(abs(del_ref)), eps);

    fprintf('\n%-28s %14s %14s %10s\n', 'quantity', 'Simulink', 'reference', 'rel diff');
    fprintf('%-28s %14.1f %14.1f %9.3f%%\n', 'peak root moment [kN m]', ...
            max(abs(mom_sl))/1e3, max(abs(mom_ref))/1e3, 100*e_mom);
    fprintf('%-28s %14.2f %14.2f %9.3f%%\n', 'peak command [deg]', ...
            rad2deg(max(abs(del_sl))), rad2deg(max(abs(del_ref))), 100*e_del);

    tol = 0.02;
    if e_mom < tol && e_del < tol
        fprintf('\nModel agrees with the reference to better than %.0f %%.\n', 100*tol);
    else
        fprintf('\nDISAGREEMENT above %.0f %%. Check Mux ordering, Selector\n', 100*tol);
        fprintf('indices and Sum signs in the model before trusting its output.\n');
    end

    figure('Name', 'Simulink against reference', 'Color', 'w');
    subplot(2,1,1);
    plot(t_sl, mom_sl/1e3, 'LineWidth', 1.6); hold on;
    plot(t_sl, mom_ref/1e3, '--', 'LineWidth', 1.4); grid on;
    ylabel('root moment [kN m]'); legend('Simulink', 'reference'); 
    title('Model verification');
    subplot(2,1,2);
    plot(t_sl, rad2deg(del_sl), 'LineWidth', 1.6); hold on;
    plot(t_sl, rad2deg(del_ref), '--', 'LineWidth', 1.4); grid on;
    xlabel('time [s]'); ylabel('command [deg]'); legend('Simulink', 'reference');
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
