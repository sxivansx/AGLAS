% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 6 : SAFETY ZONES AND THICKNESS SIZING
% =========================================================================
% Classifies the factor-of-safety field into safety zones and sizes the spar
% wall thickness to meet the target factor of safety at least mass.
%
% Corrections relative to the original version of this script:
%   - It loaded fos_data.mat expecting a variable named FoS_matrix. No stage
%     ever saved one, so the script raised "FoS_matrix undefined" and could
%     not run at all against the data in the repository.
%   - The FoS field it did plot had been overwritten one stage earlier by
%     FoS_matrix = rand(50,500)*3 + 8, commented "Dummy safe data". The
%     reported 0.270 thickness factor and 73 % mass saving were derived from
%     uniform random numbers, not from the structure.
%   - The optimiser assumed FoS scales linearly with thickness. That is wrong
%     for a solid rectangular section, where FoS goes as thickness squared,
%     and only approximate for a thin-walled box. Each candidate is now fully
%     re-analysed instead.
%   - Mass saving was reported as 1 minus the thickness factor. Mass scales
%     with cross-sectional area, not with wall thickness, so the two differ.
%   - Axes were labelled "chordwise elements" on a field whose axes are span
%     and time.
%
% Inputs  : data/strain_analysis.mat
% Outputs : data/fos_optimisation.mat, results/fos_optimisation.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'strain_analysis.mat'), ...
     'cfg', 'x', 't_sol', 'FoS_field', 'min_FoS');

fprintf('=== AGLAS supplementary stage 6: safety zones and sizing ===\n\n');

target   = cfg.safety.fos_target;
critical = cfg.safety.fos_critical;

% ------------------------------------------------- 1. safety zone map --
zones = ones(size(FoS_field));                      % 1 = critical
zones(FoS_field >= critical) = 2;                   % 2 = caution
zones(FoS_field >= target)   = 3;                   % 3 = safe
safe_mask = FoS_field >= target;

fprintf('Baseline design\n');
fprintf('  minimum FoS                 : %.4f\n', min_FoS);
fprintf('  field at or above target    : %.3f %%\n', 100*mean(safe_mask(:)));
fprintf('  field below critical        : %.3f %%\n', 100*mean(FoS_field(:) < critical));

% ------------------------------- 2. how does FoS really scale with t? --
t0 = cfg.geom.wall_t;
factors = [0.5 0.7 0.85 1.0 1.2 1.5 2.0];
scan = struct('factor', num2cell(factors));

fprintf('\nThickness scan (each point is a full re-analysis)\n');
fprintf('  %8s %10s %10s %10s %10s %9s\n', ...
        'factor', 'wall [mm]', 'min FoS', 'mass [kg]', 'f1 [Hz]', 'tip [m]');
for i = 1:numel(factors)
    c = cfg;
    c.geom.wall_t = t0 * factors(i);
    r = evaluate_design(c);
    scan(i).result = r;
    fprintf('  %8.2f %10.3f %10.4f %10.2f %10.3f %9.4f\n', ...
            factors(i), 1000*c.geom.wall_t, r.min_FoS, r.mass, r.f1, r.tip_defl);
end

fs   = [scan.factor];
fos  = arrayfun(@(s) s.result.min_FoS, scan);
mass = arrayfun(@(s) s.result.mass,    scan);
p_fit = polyfit(log(fs), log(fos), 1);
fprintf('\n  Fitted exponent, FoS proportional to thickness^n : n = %.3f\n', p_fit(1));
fprintf('  The original optimiser assumed n = 1 exactly. For this thin-walled\n');
fprintf('  box the true exponent is close to but not equal to one, and for the\n');
fprintf('  legacy solid section it would be near two.\n');

% ------------------------------------- 3. size the wall for the target --
fprintf('\nSizing the wall thickness for a minimum FoS of %.2f\n', target);

lo = 0.2; hi = 3.0;
f_lo = evaluate_design(setfield(cfg, 'geom', setfield(cfg.geom, 'wall_t', t0*lo))).min_FoS;
f_hi = evaluate_design(setfield(cfg, 'geom', setfield(cfg.geom, 'wall_t', t0*hi))).min_FoS;

if f_lo > target
    warning('AGLAS:sizingBracket', ...
        'Even the thinnest wall considered exceeds the target; no sizing needed.');
    best_factor = lo;
elseif f_hi < target
    warning('AGLAS:sizingBracket', ...
        'The target is not reachable within the thickness range considered.');
    best_factor = hi;
else
    for it = 1:30
        mid  = 0.5*(lo + hi);
        f_mid = evaluate_design(setfield(cfg, 'geom', ...
                    setfield(cfg.geom, 'wall_t', t0*mid))).min_FoS;
        if f_mid < target
            lo = mid;
        else
            hi = mid;
        end
        if (hi - lo) < 1e-4
            break;
        end
    end
    best_factor = hi;                  % the side that satisfies the target
end

cfg_opt = cfg;
cfg_opt.geom.wall_t = t0 * best_factor;
opt  = evaluate_design(cfg_opt);
base = evaluate_design(cfg);

fprintf('  thickness factor            : %.4f\n', best_factor);
fprintf('  wall thickness              : %.3f mm  ->  %.3f mm\n', 1000*t0, 1000*cfg_opt.geom.wall_t);
fprintf('  minimum FoS                 : %.4f  ->  %.4f\n', base.min_FoS, opt.min_FoS);
fprintf('  peak stress                 : %.2f MPa  ->  %.2f MPa\n', ...
        base.peak_stress/1e6, opt.peak_stress/1e6);
fprintf('  section area                : %.6f m^2  ->  %.6f m^2\n', base.area, opt.area);
fprintf('  structural mass (semi-span) : %.2f kg  ->  %.2f kg\n', base.mass, opt.mass);
fprintf('  mass change                 : %+.2f kg  (%+.2f %%)\n', ...
        opt.mass - base.mass, 100*(opt.mass - base.mass)/base.mass);
fprintf('  fundamental frequency       : %.3f Hz  ->  %.3f Hz\n', base.f1, opt.f1);
fprintf('\n  Mass follows cross-sectional area, not wall thickness. Here a %+.1f %%\n', ...
        100*(best_factor-1));
fprintf('  thickness change gives a %+.1f %% area change; the two are close for a\n', ...
        100*(opt.area-base.area)/base.area);
fprintf('  thin-walled box but not equal, and they diverge for a solid section.\n');
fprintf('  The fundamental frequency is almost unchanged because thickening the\n');
fprintf('  walls adds stiffness and mass in nearly the same proportion.\n');

% ------------------------------------------------------------- plots --
fig = figure('Name', 'Safety zones and sizing', 'NumberTitle', 'off', ...
             'Position', [100 100 1000 780]);

subplot(2,2,1);
imagesc(t_sol, x, zones); set(gca, 'YDir', 'normal');
colormap(gca, [0.85 0.20 0.20; 0.95 0.85 0.25; 0.25 0.70 0.35]);
cb = colorbar('Ticks', [1.33 2 2.66], 'TickLabels', ...
     {sprintf('< %.1f', critical), sprintf('%.1f-%.1f', critical, target), ...
      sprintf('>= %.1f', target)});
xlabel('time [s]'); ylabel('spanwise station [m]');
title('Safety zones over span and time');

subplot(2,2,2);
plot(fs, fos, '-o', 'LineWidth', 1.7); hold on;
plot(fs, fos(fs==1)*fs, '--', 'LineWidth', 1.3);
plot([fs(1) fs(end)], target*[1 1], 'k--', 'LineWidth', 1.1);
grid on; xlabel('wall thickness factor'); ylabel('minimum factor of safety');
title(sprintf('FoS vs thickness (fitted exponent %.2f)', p_fit(1)));
legend('re-analysed', 'linear assumption', sprintf('target %.1f', target), ...
       'Location', 'northwest');

subplot(2,2,3);
plot(fs, mass, '-o', 'LineWidth', 1.7); grid on;
xlabel('wall thickness factor'); ylabel('structural mass [kg]');
title('Semi-span structural mass vs thickness');

subplot(2,2,4);
plot(base.x, base.FoS_span, 'LineWidth', 1.7); hold on;
plot(opt.x,  opt.FoS_span,  'LineWidth', 1.7);
plot(base.x([1 end]), target*[1 1],   'k--', 'LineWidth', 1.1);
plot(base.x([1 end]), critical*[1 1], 'k:',  'LineWidth', 1.2);
grid on; xlabel('spanwise station [m]'); ylabel('minimum factor of safety');
title('Worst-case FoS along the span, before and after sizing');
legend('baseline', 'sized', sprintf('target %.1f', target), ...
       sprintf('critical %.1f', critical), 'Location', 'northwest');
ylim([0, cfg.safety.fos_plot_cap]);

save_figure(fig, fullfile(paths.results, 'fos_optimisation.png'), 150);

out_file = fullfile(paths.data, 'fos_optimisation.mat');
save(out_file, 'cfg', 'cfg_opt', 'zones', 'safe_mask', 'scan', 'fs', 'fos', ...
     'mass', 'p_fit', 'best_factor', 'base', 'opt', 'target', 'critical');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'fos_optimisation.png'));
