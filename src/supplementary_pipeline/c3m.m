% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 3 : GUST RESPONSE COMPARISON
% =========================================================================
% Compares the wing response under three load models:
%   1  quasi-steady aerodynamics, discrete 1-cosine gust
%   2  Kussner unsteady lag,      discrete 1-cosine gust
%   3  Kussner unsteady lag,      von Karman continuous turbulence
%
% Quantifies the dynamic amplification factor, the ratio of the peak dynamic
% deflection to the deflection the same peak load would produce statically.
%
% Inputs  : data/structural_model.mat
% Outputs : data/turbulence_response.mat, results/response_comparison.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'structural_model.mat'), 'cfg', 'fem', 'modes');

fprintf('=== AGLAS supplementary stage 3: response comparison ===\n\n');

n_modes = cfg.fem.n_modes;
md = modes;
md.Phi     = modes.Phi(:, 1:n_modes);
md.omega   = modes.omega(1:n_modes);
md.freq_hz = modes.freq_hz(1:n_modes);
md.n_modes = n_modes;

% --- Case 1 and 2: discrete gust ---------------------------------------
t1 = 0:cfg.time.dt_out:cfg.time.t_end;
wg1 = gust_one_minus_cos(t1, cfg.gust.U_ds, cfg.gust.t_g);

[tA, qA, UA, giA] = solve_modal_response(cfg, fem, md, t1, wg1, 'quasisteady', cfg.gust.t_g/20);
[tB, qB, UB, giB] = solve_modal_response(cfg, fem, md, t1, wg1, 'kussner',     cfg.gust.t_g/20);
tipA = UA(:, fem.idx_w(end));
tipB = UB(:, fem.idx_w(end));

% Static deflection under the same peak load, for the amplification factor.
f_static  = consistent_line_load(fem, max(abs(giB.q_line)), max(abs(giB.m_line)));
u_static  = fem.K \ f_static;
w_static  = u_static(fem.idx_w(end));

fprintf('Discrete 1-cosine gust, %0.1f m/s over %.1f s\n', cfg.gust.U_ds, cfg.gust.t_g);
fprintf('  peak tip deflection, quasi-steady : %8.4f m\n', max(abs(tipA)));
fprintf('  peak tip deflection, Kussner      : %8.4f m\n', max(abs(tipB)));
fprintf('  Kussner reduces the peak by        %8.2f %%\n', ...
        100*(max(abs(tipA))-max(abs(tipB)))/max(abs(tipA)));
fprintf('  static deflection at peak load    : %8.4f m\n', w_static);
fprintf('  dynamic amplification factor      : %8.3f\n', max(abs(tipB))/w_static);
fprintf('  first bending period %.3f s vs gust duration %.2f s\n', ...
        1/md.freq_hz(1), cfg.gust.t_g);
fprintf('  The gust is slow relative to the structure, so the amplification\n');
fprintf('  factor is close to one, as expected.\n');

% --- Case 3: continuous turbulence -------------------------------------
T_turb = 30;   % long enough for stable RMS statistics, short enough to run fast
[wg3, t3] = gust_von_karman(cfg, T_turb, cfg.time.dt_out);
[tC, qC, UC, giC] = solve_modal_response(cfg, fem, md, t3, wg3, 'kussner', 0.02);
tipC = UC(:, fem.idx_w(end));

fprintf('\nvon Karman continuous turbulence, %.0f s record\n', T_turb);
fprintf('  gust RMS                          : %8.3f m/s\n', std(wg3));
fprintf('  tip deflection RMS                : %8.4f m\n', std(tipC));
fprintf('  tip deflection peak               : %8.4f m\n', max(abs(tipC)));
fprintf('  peak-to-RMS ratio                 : %8.2f\n', max(abs(tipC))/std(tipC));

% --- Plot ---------------------------------------------------------------
fig = figure('Name', 'Gust response comparison', 'NumberTitle', 'off', ...
             'Position', [100 100 950 900]);

subplot(4,1,1);
plot(t1, giA.q_line/1000, 'LineWidth', 1.6); hold on;
plot(t1, giB.q_line/1000, '--', 'LineWidth', 1.6);
grid on; ylabel('line load [kN/m]');
legend('quasi-steady', 'Kussner lag', 'Location', 'northeast');
title('Aerodynamic line load from the discrete gust');

subplot(4,1,2);
plot(tA, tipA, 'LineWidth', 1.6); hold on;
plot(tB, tipB, '--', 'LineWidth', 1.6);
plot([tA(1) tA(end)], [w_static w_static], 'k:', 'LineWidth', 1.2);
grid on; xlabel('time [s]'); ylabel('tip deflection [m]');
legend('quasi-steady', 'Kussner lag', 'static at peak load', 'Location', 'northeast');
title('Tip response to the discrete gust');

subplot(4,1,3);
plot(t3, wg3, 'LineWidth', 0.9); grid on;
ylabel('gust [m/s]'); xlim([t3(1) t3(end)]);
title(sprintf('von Karman turbulence input (RMS %.2f m/s)', std(wg3)));

subplot(4,1,4);
plot(tC, tipC, 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('tip deflection [m]'); xlim([t3(1) t3(end)]);
title(sprintf('Tip response to turbulence (RMS %.3f m, peak %.3f m)', ...
              std(tipC), max(abs(tipC))));

save_figure(fig, fullfile(paths.results, 'response_comparison.png'), 150);

out_file = fullfile(paths.data, 'turbulence_response.mat');
save(out_file, 'cfg', 'tA', 'tipA', 'tB', 'tipB', 'tC', 'tipC', ...
     'qA', 'qB', 'qC', 'wg1', 'wg3', 't1', 't3', 'w_static');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'response_comparison.png'));
