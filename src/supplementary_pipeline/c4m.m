% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 4 : BENDING STRESS FIELD
% =========================================================================
% Recovers the spanwise and time-resolved bending moment and stress fields
% from the gust response, and cross-checks the peak against a static hand
% calculation.
%
% Corrections relative to the original version of this script:
%   - The bending moment was taken to be Phi*q, which is the displacement
%     field, and then divided by I to get a stress. That is dimensionally
%     invalid: it produces metres per metre to the fourth, not pascals.
%     Moment is EI times curvature, so a second spatial derivative is
%     required and is now taken exactly from the element interpolation.
%   - The mode shapes were hand-written as sin(i*pi*x/(2L)), which has a
%     non-zero slope at the clamped root and so violates the cantilever
%     boundary condition. The real FE modes are used instead.
%   - The modal response was fabricated as 0.02*sin(2*pi*i*t)*exp(-0.3*t)
%     rather than taken from the simulation, so nothing downstream depended
%     on the actual physics.
%   - The script contradicted itself on geometry, using I = 5e-6 and then
%     I = b*h^3/12 = 1e-7, and z = 0.05 m for a section only 0.02 m deep.
%     All geometry now comes from aglas_config.
%
% Inputs  : data/structural_model.mat, data/modal_response.mat
% Outputs : data/stress_field.mat, results/stress_field.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'structural_model.mat'), 'fem');
load(fullfile(paths.data, 'modal_response.mat'), 'cfg', 't_sol', 'U_dof', 'gust_info');

fprintf('=== AGLAS supplementary stage 4: bending stress field ===\n\n');

rec = recover_bending(cfg, fem, U_dof.');

x       = rec.x(:);
moment  = rec.moment;              % [n_nodes x n_time]
stress  = rec.stress;
sigma_y = cfg.material.sigma_yield;

[peak_M, iM]  = max(abs(moment(:)));
[nM, tM]      = ind2sub(size(moment), iM);
[peak_s, is]  = max(abs(stress(:)));
[ns, ts]      = ind2sub(size(stress), is);

fprintf('Peak bending moment : %10.1f N.m at x = %.2f m, t = %.3f s\n', ...
        peak_M, x(nM), t_sol(tM));
fprintf('Peak bending stress : %10.2f MPa at x = %.2f m, t = %.3f s\n', ...
        peak_s/1e6, x(ns), t_sol(ts));
fprintf('Yield strength      : %10.2f MPa\n', sigma_y/1e6);
fprintf('Stress utilisation  : %10.1f %% of yield\n', 100*peak_s/sigma_y);

% --- Independent static cross-check -------------------------------------
q_peak = max(abs(gust_info.q_line));
M_hand = q_peak * cfg.geom.L^2 / 2;
s_hand = M_hand * cfg.section.c_outer / cfg.section.I;

fprintf('\nStatic cross-check at the peak line load of %.1f N/m\n', q_peak);
fprintf('  hand calculation  M = q*L^2/2      : %10.1f N.m\n', M_hand);
fprintf('  hand calculation  sigma = M*c/I    : %10.2f MPa\n', s_hand/1e6);
fprintf('  dynamic result relative to static  : %10.3f\n', peak_s/s_hand);
fprintf('  The ratio is the dynamic amplification factor. A value near one is\n');
fprintf('  expected because the %.1f s gust is slow next to the %.3f s period\n', ...
        cfg.gust.t_g, 1/sqrt(cfg.section.EI/(cfg.section.m_total*cfg.geom.L^4))/0.5596/(2*pi)*(2*pi));

% --- Plot ---------------------------------------------------------------
fig = figure('Name', 'Bending stress field', 'NumberTitle', 'off', ...
             'Position', [100 100 1000 760]);

subplot(2,2,1);
plot(x, moment(:, tM)/1000, 'LineWidth', 1.8); hold on;
plot(x, (q_peak*(cfg.geom.L - x).^2/2)/1000, '--', 'LineWidth', 1.4);
grid on; xlabel('spanwise station [m]'); ylabel('bending moment [kN m]');
title('Spanwise bending moment at the peak instant');
legend('FE recovery', 'static q(L-x)^2/2', 'Location', 'northeast');

subplot(2,2,2);
plot(t_sol, stress(1,:)/1e6, 'LineWidth', 1.6); hold on;
plot(t_sol, stress(round(end/2),:)/1e6, '--', 'LineWidth', 1.4);
grid on; xlabel('time [s]'); ylabel('stress [MPa]');
title('Stress history');
legend('root', 'mid-span', 'Location', 'northeast');

subplot(2,2,3);
imagesc(t_sol, x, stress/1e6); set(gca, 'YDir', 'normal');
colorbar; xlabel('time [s]'); ylabel('spanwise station [m]');
title('Stress field [MPa]');

subplot(2,2,4);
surf(t_sol, x, stress/1e6, 'EdgeColor', 'none');
xlabel('time [s]'); ylabel('span [m]'); zlabel('stress [MPa]');
title('Stress over span and time'); view(40, 32); colorbar;

save_figure(fig, fullfile(paths.results, 'stress_field.png'), 150);

out_file = fullfile(paths.data, 'stress_field.mat');
save(out_file, 'cfg', 'x', 't_sol', 'moment', 'stress', 'peak_s', 'peak_M', ...
     'M_hand', 's_hand');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'stress_field.png'));
