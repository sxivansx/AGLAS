% =========================================================================
% AGLAS STAGE 3 (main pipeline) : GUST RESPONSE
% =========================================================================
% Integrates the modal equations of motion under a discrete "1-cosine" gust
% and reconstructs the physical deflection field.
%
% Modal equations, with mass-normalised mode shapes:
%     qddot_i + 2*zeta_i*omega_i*qdot_i + omega_i^2*q_i = F_i(t)
%
% The generalised force F(t) comes from strip-theory aerodynamics with a
% Kussner gust-penetration lag. See gust_modal_force for the derivation and
% for what the original implementation got wrong.
%
% Inputs  : data/structural_model.mat
% Outputs : data/modal_response.mat, results/gust_response.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'structural_model.mat'), 'cfg', 'fem', 'modes');

fprintf('=== AGLAS Stage 3: gust response ===\n\n');

n_modes = cfg.fem.n_modes;
md = modes;
md.Phi     = modes.Phi(:, 1:n_modes);
md.omega   = modes.omega(1:n_modes);
md.freq_hz = modes.freq_hz(1:n_modes);
md.n_modes = n_modes;

fprintf('Retaining %d modes: %s Hz\n', n_modes, ...
        mat2str(round(md.freq_hz(:).'*1000)/1000));

% --- Gust and generalised force ----------------------------------------
t_grid = 0:cfg.time.dt_out:cfg.time.t_end;
wg     = gust_one_minus_cos(t_grid, cfg.gust.U_ds, cfg.gust.t_g);

[t_sol, q_sol, U_dof, gust_info] = ...
    solve_modal_response(cfg, fem, md, t_grid, wg, 'kussner', cfg.gust.t_g/20);

fprintf('Gust: 1-cosine, peak %.1f m/s over %.2f s at U = %.1f m/s\n', ...
        cfg.gust.U_ds, cfg.gust.t_g, cfg.flight.U_inf);
fprintf('Peak incidence change  : %.2f deg\n', gust_info.d_alpha_deg);
fprintf('Peak lift line load    : %.1f N/m  (%.0f N on the semi-span)\n', ...
        max(abs(gust_info.q_line)), gust_info.peak_lift_N);
if ~gust_info.linear_aero_ok
    warning('AGLAS:linearAero', ...
      ['Peak incidence change is %.2f deg, above the %.1f deg linear-aerodynamics ' ...
       'limit. Lift is assumed proportional to incidence, so the computed loads ' ...
       'are an upper bound; a real aerofoil would begin to stall.'], ...
       gust_info.d_alpha_deg, cfg.limits.max_delta_alpha);
end

% --- Reconstruct the physical field -------------------------------------
w_field = [zeros(size(q_sol,1),1), U_dof(:, fem.idx_w)];
tip_w   = w_field(:, end);

[peak_tip, i_peak] = max(abs(tip_w));
tip_frac = peak_tip / cfg.geom.L;

fprintf('\nPeak tip deflection    : %.4f m at t = %.3f s  (%.2f %% of span)\n', ...
        peak_tip, t_sol(i_peak), 100*tip_frac);
if tip_frac > cfg.limits.max_tip_defl_frac
    warning('AGLAS:largeDeflection', ...
      ['Tip deflection is %.1f %% of span, beyond the %.0f %% small-deflection ' ...
       'limit of linear beam theory. Geometric stiffening is neglected, so the ' ...
       'computed deflection overstates the real one.'], ...
       100*tip_frac, 100*cfg.limits.max_tip_defl_frac);
end

if fem.torsion
    tw = [zeros(size(q_sol,1),1), U_dof(:, fem.idx_ph)];
    fprintf('Peak tip twist         : %.4f deg\n', max(abs(tw(:,end)))*180/pi);
else
    tw = [];
end

% --- Plot ---------------------------------------------------------------
fig = figure('Name', 'Gust response', 'NumberTitle', 'off', ...
             'Position', [100 100 860 760]);

subplot(3,1,1);
plot(t_grid, wg, 'LineWidth', 1.6); hold on;
plot(t_grid, gust_info.wg_effective, '--', 'LineWidth', 1.4);
grid on; ylabel('gust velocity [m/s]');
legend('1-cosine gust', 'after Kussner lag', 'Location', 'northeast');
title('Gust input');

subplot(3,1,2);
plot(t_sol, tip_w, 'LineWidth', 1.6); hold on;
plot(t_sol(i_peak), tip_w(i_peak), 'o', 'MarkerSize', 7, 'LineWidth', 1.4);
grid on; ylabel('tip deflection [m]');
title(sprintf('Wing tip response  (peak %.4f m at t = %.2f s)', ...
              tip_w(i_peak), t_sol(i_peak)));

subplot(3,1,3);
plot(t_sol, q_sol, 'LineWidth', 1.3);
grid on; xlabel('time [s]'); ylabel('modal coordinate');
lg = cell(1, n_modes);
for i = 1:n_modes
    lg{i} = sprintf('mode %d (%.2f Hz)', i, md.freq_hz(i));
end
legend(lg, 'Location', 'northeast');
title('Generalised coordinates');

save_figure(fig, fullfile(paths.results, 'gust_response.png'), 150);

out_file = fullfile(paths.data, 'modal_response.mat');
save(out_file, 'cfg', 't_sol', 'q_sol', 'U_dof', 'w_field', 'tip_w', ...
     'tw', 'n_modes', 'gust_info', 'wg', 't_grid');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'gust_response.png'));
