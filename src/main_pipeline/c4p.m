% =========================================================================
% AGLAS STAGE 4 (main pipeline) : FLUTTER AND DIVERGENCE
% =========================================================================
% Aeroelastic stability by the p-k method with Theodorsen strip-theory
% unsteady aerodynamics, producing the classical V-g and V-f diagrams.
%
% This analysis requires a torsional degree of freedom. A bending-only beam
% has no mechanism for the bending and torsion motions to exchange energy
% through the airstream, so it cannot flutter at any speed. The original
% version of this script ran on a bending-only model with a hand-written,
% real, symmetric positive-definite "aerodynamic" matrix, which could only
% ever add damping; every mode stayed stable at every speed and the script
% saved flutter_speed = NaN. The torsion DOF added in beam_fem is what makes
% a flutter speed exist at all.
%
% Inputs  : data/structural_model.mat
% Outputs : data/aero_data.mat, results/flutter_vgf.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'structural_model.mat'), 'cfg', 'fem', 'modes');

fprintf('=== AGLAS Stage 4: flutter and divergence ===\n\n');

if ~fem.torsion
    error('AGLAS:noTorsion', ...
      ['This model has no torsional DOF, so bending-torsion flutter cannot ' ...
       'occur and a flutter speed is not defined. Set ' ...
       'cfg.fem.include_torsion = true in aglas_config.']);
end

n_fl = min(cfg.flutter.n_modes, modes.n_modes);
md = modes;
md.Phi     = modes.Phi(:, 1:n_fl);
md.omega   = modes.omega(1:n_fl);
md.freq_hz = modes.freq_hz(1:n_fl);
md.type    = modes.type(1:n_fl);
md.n_modes = n_fl;

fprintf('Modes carried into the p-k solution (%d):\n', n_fl);
for i = 1:n_fl
    fprintf('  %2d : %9.3f Hz  %s\n', i, md.freq_hz(i), md.type{i});
end

fprintf('\nSweeping %.0f to %.0f m/s in %.0f m/s steps ...\n', ...
        cfg.flutter.U(1), cfg.flutter.U(end), cfg.flutter.U(2)-cfg.flutter.U(1));

res = flutter_pk(cfg, fem, md);

fprintf('p-k iteration converged at %.1f %% of sweep points\n', ...
        100*mean(res.converged(:)));

% ------------------------------------------------------------- results --
fprintf('\n--- Stability summary ---\n');
if isnan(res.flutter_speed)
    fprintf('Flutter      : none found below %.0f m/s\n', cfg.flutter.U(end));
else
    fprintf('Flutter      : %8.2f m/s   %7.3f Hz   (branch %d)\n', ...
            res.flutter_speed, res.flutter_freq_hz, res.flutter_mode);
end
if isnan(res.divergence_speed)
    fprintf('Divergence   : none found\n');
else
    fprintf('Divergence   : %8.2f m/s\n', res.divergence_speed);
end

fprintf(['\nBranch labels are assigned by matching eigenvector shape and\n' ...
         'eigenvalue continuity between speeds. Where two branches coalesce,\n' ...
         'which is what produces classical bending-torsion flutter, their\n' ...
         'identities genuinely merge and the labels may exchange beyond that\n' ...
         'point. The flutter speed above does not depend on the labelling: it\n' ...
         'is taken from the largest real part over all roots at each speed.\n']);

U_crit = min([res.flutter_speed, res.divergence_speed]);
if isfinite(U_crit)
    fprintf('Critical speed: %.2f m/s, Mach %.3f at ISA sea level\n', ...
            U_crit, U_crit/cfg.flight.a_sound);
    fprintf('Margin over the %.0f m/s reference cruise: %.2f\n', ...
            cfg.flight.U_inf, U_crit/cfg.flight.U_inf);
    if U_crit/cfg.flight.a_sound > 0.3
        warning('AGLAS:compressibility', ...
          ['The critical speed is Mach %.2f. Theodorsen theory is ' ...
           'incompressible and is only reliable below about Mach 0.3. Treat ' ...
           'this speed as an indicative result; a compressible unsteady ' ...
           'method, or at minimum a Prandtl-Glauert correction, is needed ' ...
           'to quote it as a certified number.'], U_crit/cfg.flight.a_sound);
    end
end

% ---------------------------------------------------------------- plots --
fig = figure('Name', 'Flutter V-g-f diagram', 'NumberTitle', 'off', ...
             'Position', [100 100 880 760]);

subplot(2,1,1);
h_f = plot(res.U, res.freq_hz, 'LineWidth', 1.5);
grid on; ylabel('frequency [Hz]');
title('V-f diagram: aeroelastic frequency vs airspeed');
hold on;
if isfinite(res.flutter_speed)
    yl = ylim;
    plot([res.flutter_speed res.flutter_speed], yl, 'r--', 'LineWidth', 1.5);
    ylim(yl);
end
% Pass the branch handles explicitly so the flutter marker stays out of the
% legend. Excluding it via the Annotation property is MATLAB-only.
lg = cell(1, numel(h_f));
for i = 1:numel(lg)
    lg{i} = sprintf('branch %d (%.1f Hz)', i, res.freq_hz(1,i));
end
legend(h_f, lg, 'Location', 'northwest');

subplot(2,1,2);
plot(res.U, res.g, 'LineWidth', 1.5); hold on;
plot(res.U([1 end]), [0 0], 'k-', 'LineWidth', 1.0);
grid on; xlabel('airspeed [m/s]'); ylabel('damping g = 2 Re(p)/|Im(p)|');
title('V-g diagram: positive damping means flutter');
if isfinite(res.flutter_speed)
    yl = ylim;
    plot([res.flutter_speed res.flutter_speed], yl, 'r--', 'LineWidth', 1.5);
    text(res.flutter_speed, yl(2), ...
         sprintf('  flutter %.1f m/s, %.2f Hz', res.flutter_speed, res.flutter_freq_hz), ...
         'VerticalAlignment', 'top', 'Color', 'r');
    ylim(yl);
end
gvis = res.g(isfinite(res.g));
if ~isempty(gvis)
    ylim([max(min(gvis), -1.0), min(max(max(gvis), 0.05), 0.5)]);
end

save_figure(fig, fullfile(paths.results, 'flutter_vgf.png'), 150);

U               = res.U;
damping         = res.damping;
frequency       = res.freq_hz;
flutter_speed   = res.flutter_speed;
flutter_freq    = res.flutter_freq_hz;
divergence_speed= res.divergence_speed;

out_file = fullfile(paths.data, 'aero_data.mat');
save(out_file, 'cfg', 'res', 'U', 'damping', 'frequency', ...
     'flutter_speed', 'flutter_freq', 'divergence_speed');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'flutter_vgf.png'));
