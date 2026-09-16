% =========================================================================
% AGLAS STAGE 6 (main pipeline) : STRESS, STRAIN AND FACTOR OF SAFETY
% =========================================================================
% Recovers the bending moment, outer-fibre stress and strain fields over the
% whole span and the whole gust event, then evaluates the factor of safety.
%
% Curvature comes from the exact second derivative of the cubic Hermite
% element interpolation, using both the translation and the rotation degrees
% of freedom. See recover_bending for why the original central-difference
% approach missed the root, where a cantilever carries its largest moment,
% and so reported a factor of safety 17 % higher than the true value.
%
% Inputs  : data/modal_response.mat, data/structural_model.mat
% Outputs : data/fos_data.mat, results/stress_analysis.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'structural_model.mat'), 'fem');
load(fullfile(paths.data, 'modal_response.mat'), 'cfg', 't_sol', 'U_dof');

fprintf('=== AGLAS Stage 6: stress and factor of safety ===\n\n');

rec = recover_bending(cfg, fem, U_dof.');     % [n_nodes x n_time]

x          = rec.x(:);
stress     = rec.stress;
strain     = rec.strain;
moment     = rec.moment;
sigma_y    = cfg.material.sigma_yield;

FoS_field  = sigma_y ./ max(abs(stress), eps);

[max_stress, lin_idx]       = max(abs(stress(:)));
[i_node_max, i_time_max]    = ind2sub(size(stress), lin_idx);
min_FoS                     = sigma_y / max_stress;

root_stress_t               = stress(1, :);
[peak_root_stress, i_root_t]= max(abs(root_stress_t));
FoS_root                    = sigma_y / peak_root_stress;

fprintf('Peak bending moment    : %10.1f N.m  at x = %.2f m\n', ...
        max(abs(moment(:))), x(i_node_max));
fprintf('Peak stress            : %10.2f MPa at x = %.2f m, t = %.3f s\n', ...
        max_stress/1e6, x(i_node_max), t_sol(i_time_max));
fprintf('Peak strain            : %10.2f microstrain\n', max(abs(strain(:)))*1e6);
fprintf('Peak stress at the root: %10.2f MPa at t = %.3f s\n', ...
        peak_root_stress/1e6, t_sol(i_root_t));
fprintf('Yield strength         : %10.2f MPa (%s)\n', sigma_y/1e6, cfg.material.name);
fprintf('\nMinimum factor of safety : %.3f  (at x = %.2f m)\n', min_FoS, x(i_node_max));
fprintf('Factor of safety at root : %.3f\n', FoS_root);

fprintf('\nAcceptance\n');
if min_FoS < cfg.safety.fos_critical
    fprintf('  FAIL  : below the critical threshold of %.1f\n', cfg.safety.fos_critical);
elseif min_FoS < cfg.safety.fos_target
    fprintf('  MARGINAL : above the critical %.1f but below the target %.1f\n', ...
            cfg.safety.fos_critical, cfg.safety.fos_target);
else
    fprintf('  PASS  : at or above the target of %.1f\n', cfg.safety.fos_target);
end

fprintf('\nMesh quality: peak inter-element curvature jump %.3e 1/m\n', ...
        max(abs(rec.jump(:))));

% ---------------------------------------------------------------- plots --
fig = figure('Name', 'Stress and factor of safety', 'NumberTitle', 'off', ...
             'Position', [100 100 1000 780]);

subplot(2,2,1);
plot(t_sol, root_stress_t/1e6, 'LineWidth', 1.6); hold on;
plot(t_sol(i_root_t), root_stress_t(i_root_t)/1e6, 'o', 'MarkerSize', 7, 'LineWidth', 1.4);
grid on; xlabel('time [s]'); ylabel('stress [MPa]');
title(sprintf('Root stress history (peak %.1f MPa)', peak_root_stress/1e6));

subplot(2,2,2);
plot(x, stress(:, i_time_max)/1e6, 'LineWidth', 1.8); hold on;
plot(x([1 end]), [ sigma_y  sigma_y]/1e6, 'k--', 'LineWidth', 1.0);
plot(x([1 end]), [-sigma_y -sigma_y]/1e6, 'k--', 'LineWidth', 1.0);
grid on; xlabel('spanwise station [m]'); ylabel('stress [MPa]');
title(sprintf('Spanwise stress at peak, t = %.2f s', t_sol(i_time_max)));
legend('stress', 'yield', 'Location', 'northeast');

subplot(2,2,3);
imagesc(t_sol, x, min(FoS_field, cfg.safety.fos_plot_cap));
set(gca, 'YDir', 'normal');
colorbar; xlabel('time [s]'); ylabel('spanwise station [m]');
title(sprintf('Factor of safety field (clipped at %g)', cfg.safety.fos_plot_cap));

subplot(2,2,4);
plot(x, min(FoS_field, [], 2), 'LineWidth', 1.8); hold on;
plot(x([1 end]), cfg.safety.fos_target*[1 1],   'k--', 'LineWidth', 1.0);
plot(x([1 end]), cfg.safety.fos_critical*[1 1], 'k:',  'LineWidth', 1.2);
grid on; xlabel('spanwise station [m]'); ylabel('minimum factor of safety');
title('Minimum factor of safety along the span');
legend('minimum FoS', sprintf('target %.1f', cfg.safety.fos_target), ...
       sprintf('critical %.1f', cfg.safety.fos_critical), 'Location', 'northwest');
ylim([0, min(cfg.safety.fos_plot_cap, max(min(FoS_field, [], 2))*1.2)]);

save_figure(fig, fullfile(paths.results, 'stress_analysis.png'), 150);

out_file = fullfile(paths.data, 'fos_data.mat');
save(out_file, 'cfg', 'x', 't_sol', 'stress', 'strain', 'moment', ...
     'FoS_field', 'min_FoS', 'FoS_root', 'max_stress', 'peak_root_stress', ...
     'sigma_y');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'stress_analysis.png'));
