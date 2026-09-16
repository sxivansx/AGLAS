% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 5 : STRAIN AND SAFETY VALIDATION
% =========================================================================
% Converts the recovered stress field into strain, evaluates the factor of
% safety over span and time, and reports where and when the structure is
% most highly loaded.
%
% Corrections relative to the original version of this script:
%   - The spanwise grid was built as linspace(0, L, N) with N the number of
%     ELEMENTS, while the mode shape array held one row per NODE excluding
%     the clamped root. The grid spacing was therefore wrong by L/N against
%     L/(N+1), and the first row of the mode shapes, which sits one element
%     outboard of the root, was placed at x = 0. Curvature goes as 1/dx^2,
%     so an 11 % error in dx became a 23 % error in every stress and every
%     factor of safety the script reported.
%   - The header promised strain_analysis_results.mat, which was never
%     written.
%   - caxis was renamed clim in R2022a and WindowState 'maximized' fails in
%     headless sessions; neither is used now.
%
% Inputs  : data/stress_field.mat
% Outputs : data/strain_analysis.mat, results/strain_safety.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'stress_field.mat'), 'cfg', 'x', 't_sol', 'stress');

fprintf('=== AGLAS supplementary stage 5: strain and safety ===\n\n');

E       = cfg.material.E;
sigma_y = cfg.material.sigma_yield;

strain     = stress / E;
FoS_field  = sigma_y ./ max(abs(stress), eps);
FoS_span   = min(FoS_field, [], 2);           % worst case in time, per station

[max_strain, k]        = max(abs(strain(:)));
[i_node, i_time]       = ind2sub(size(strain), k);
[min_FoS, kf]          = min(FoS_field(:));
[jf_node, jf_time]     = ind2sub(size(FoS_field), kf);

fprintf('Peak strain          : %9.1f microstrain at x = %.2f m, t = %.3f s\n', ...
        max_strain*1e6, x(i_node), t_sol(i_time));
fprintf('Peak stress          : %9.2f MPa\n', max(abs(stress(:)))/1e6);
fprintf('Yield strain         : %9.1f microstrain\n', sigma_y/E*1e6);
fprintf('Minimum FoS          : %9.3f at x = %.2f m, t = %.3f s\n', ...
        min_FoS, x(jf_node), t_sol(jf_time));

fprintf('\nSafety zones (thresholds %.1f critical / %.1f target)\n', ...
        cfg.safety.fos_critical, cfg.safety.fos_target);
n_tot  = numel(FoS_field);
n_crit = sum(FoS_field(:) <  cfg.safety.fos_critical);
n_caut = sum(FoS_field(:) >= cfg.safety.fos_critical & FoS_field(:) < cfg.safety.fos_target);
n_safe = sum(FoS_field(:) >= cfg.safety.fos_target);
fprintf('  critical  (FoS < %.1f)      : %7.3f %% of the span-time field\n', ...
        cfg.safety.fos_critical, 100*n_crit/n_tot);
fprintf('  caution   (%.1f to %.1f)     : %7.3f %%\n', ...
        cfg.safety.fos_critical, cfg.safety.fos_target, 100*n_caut/n_tot);
fprintf('  safe      (FoS >= %.1f)     : %7.3f %%\n', ...
        cfg.safety.fos_target, 100*n_safe/n_tot);

fprintf('\nVerdict: ');
if min_FoS < cfg.safety.fos_critical
    fprintf('FAIL, minimum FoS %.3f is below the critical threshold.\n', min_FoS);
elseif min_FoS < cfg.safety.fos_target
    fprintf('MARGINAL, minimum FoS %.3f clears the critical threshold but not\n', min_FoS);
    fprintf('         the %.1f target. Stage 6 sizes the section to close the gap.\n', ...
            cfg.safety.fos_target);
else
    fprintf('PASS, minimum FoS %.3f meets the %.1f target.\n', min_FoS, cfg.safety.fos_target);
end

% ------------------------------------------------------------- plots --
fig = figure('Name', 'Strain and safety', 'NumberTitle', 'off', ...
             'Position', [100 100 1000 780]);

subplot(2,2,1);
plot(t_sol, strain(1,:)*1e6, 'LineWidth', 1.6); grid on;
xlabel('time [s]'); ylabel('strain [microstrain]');
title('Root strain history');

subplot(2,2,2);
plot(x, strain(:, i_time)*1e6, 'LineWidth', 1.8); grid on;
xlabel('spanwise station [m]'); ylabel('strain [microstrain]');
title(sprintf('Spanwise strain at t = %.2f s', t_sol(i_time)));

subplot(2,2,3);
imagesc(t_sol, x, min(FoS_field, cfg.safety.fos_plot_cap));
set(gca, 'YDir', 'normal'); colorbar;
xlabel('time [s]'); ylabel('spanwise station [m]');
title(sprintf('Factor of safety field (clipped at %g)', cfg.safety.fos_plot_cap));
hold on; plot(t_sol(jf_time), x(jf_node), 'w+', 'MarkerSize', 12, 'LineWidth', 2);

subplot(2,2,4);
plot(x, FoS_span, 'LineWidth', 1.8); hold on;
plot(x([1 end]), cfg.safety.fos_target*[1 1],   'k--', 'LineWidth', 1.1);
plot(x([1 end]), cfg.safety.fos_critical*[1 1], 'k:',  'LineWidth', 1.3);
grid on; xlabel('spanwise station [m]'); ylabel('minimum factor of safety');
title('Worst-case factor of safety along the span');
legend('minimum FoS', sprintf('target %.1f', cfg.safety.fos_target), ...
       sprintf('critical %.1f', cfg.safety.fos_critical), 'Location', 'northwest');
ylim([0, cfg.safety.fos_plot_cap]);

save_figure(fig, fullfile(paths.results, 'strain_safety.png'), 150);

out_file = fullfile(paths.data, 'strain_analysis.mat');
save(out_file, 'cfg', 'x', 't_sol', 'strain', 'stress', 'FoS_field', ...
     'FoS_span', 'min_FoS', 'max_strain');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'strain_safety.png'));
