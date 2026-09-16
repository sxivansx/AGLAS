% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 1 : GUST MODELS
% =========================================================================
% Generates and compares the two atmospheric disturbance models used by the
% project: the CS-25.341 discrete "1-cosine" gust and von Karman continuous
% turbulence. Verifies the turbulence realisation against its target
% spectrum.
%
% Outputs : data/gust_models.mat, results/gust_models.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

cfg = aglas_config('baseline');

fprintf('=== AGLAS supplementary stage 1: gust models ===\n\n');

% --- Discrete 1-cosine gust --------------------------------------------
dt = cfg.time.dt_out;
t  = 0:dt:cfg.time.t_end;
wg_disc = gust_one_minus_cos(t, cfg.gust.U_ds, cfg.gust.t_g);

fprintf('Discrete 1-cosine gust\n');
fprintf('  peak velocity        : %.2f m/s at t = %.3f s\n', ...
        max(wg_disc), t(find(wg_disc == max(wg_disc), 1)));
fprintf('  duration             : %.2f s\n', cfg.gust.t_g);
fprintf('  gust gradient length : %.1f m at U = %.1f m/s\n', ...
        cfg.flight.U_inf*cfg.gust.t_g/2, cfg.flight.U_inf);
fprintf('  value at t = 0       : %.3e  (must be 0)\n', wg_disc(1));
fprintf('  value at t = t_g     : %.3e  (must be 0)\n', ...
        interp1(t, wg_disc, cfg.gust.t_g));

% --- von Karman continuous turbulence ----------------------------------
T_turb = 200;                        % long record for a converged spectrum
[wg_turb, t_turb, vk] = gust_von_karman(cfg, T_turb, dt);

fprintf('\nvon Karman continuous turbulence\n');
fprintf('  target RMS           : %.3f m/s\n', vk.rms_target);
fprintf('  realised RMS         : %.3f m/s  (%.1f %% of target)\n', ...
        vk.rms_actual, 100*vk.rms_actual/vk.rms_target);
fprintf('  mean                 : %.3e m/s  (must be ~0)\n', mean(wg_turb));
fprintf('  length scale         : %.0f m\n', cfg.turb.L_scale);
fprintf('  record length        : %.0f s at U = %.1f m/s (%.0f m of air)\n', ...
        T_turb, cfg.flight.U_inf, T_turb*cfg.flight.U_inf);
fprintf('  RNG seed             : %d (reproducible)\n', cfg.turb.seed);

% The realised RMS falls slightly short of the target because a finite record
% cannot resolve frequencies below 2*pi/T, and the von Karman spectrum still
% carries energy there. Lengthening the record closes the gap.

% --- Verify the inertial-subrange slope ---------------------------------
mask = vk.omega > 5 & vk.omega < 60;
pfit = polyfit(log(vk.omega(mask)), log(vk.S(mask)), 1);
fprintf('  PSD slope (inertial) : %.4f  (von Karman theory: -5/3 = -1.6667)\n', pfit(1));

% --- Plot ---------------------------------------------------------------
fig = figure('Name', 'Gust models', 'NumberTitle', 'off', ...
             'Position', [100 100 950 760]);

subplot(3,1,1);
plot(t, wg_disc, 'LineWidth', 1.8); grid on;
xlabel('time [s]'); ylabel('gust velocity [m/s]');
title(sprintf('Discrete 1-cosine gust: %.0f m/s peak over %.1f s', ...
              cfg.gust.U_ds, cfg.gust.t_g));

subplot(3,1,2);
plot(t_turb(t_turb <= 40), wg_turb(t_turb <= 40), 'LineWidth', 1.0); grid on;
xlabel('time [s]'); ylabel('gust velocity [m/s]');
title(sprintf('von Karman turbulence: RMS %.2f m/s, L = %.0f m (first 40 s)', ...
              vk.rms_actual, cfg.turb.L_scale));

subplot(3,1,3);
loglog(vk.omega, vk.S, 'LineWidth', 1.6); hold on;
ref = vk.S(find(mask, 1)) * (vk.omega(mask)/vk.omega(find(mask,1))).^(-5/3);
loglog(vk.omega(mask), ref, '--', 'LineWidth', 1.4);
grid on; xlabel('\omega [rad/s]'); ylabel('S_w(\omega) [(m/s)^2/(rad/s)]');
title('von Karman vertical velocity spectrum');
legend('target PSD', '\omega^{-5/3} reference', 'Location', 'southwest');

save_figure(fig, fullfile(paths.results, 'gust_models.png'), 150);

out_file = fullfile(paths.data, 'gust_models.mat');
save(out_file, 'cfg', 't', 'wg_disc', 't_turb', 'wg_turb', 'vk');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'gust_models.png'));
