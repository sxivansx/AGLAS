% =========================================================================
% AGLAS STAGE 2 (main pipeline) : STRUCTURAL MODEL AND MODAL ANALYSIS
% =========================================================================
% Assembles the cantilever finite element model, solves the free-vibration
% eigenproblem and validates the result against the closed-form
% Euler-Bernoulli cantilever solution.
%
% Inputs  : data/wing_geom.mat
% Outputs : data/structural_model.mat, results/mode_shapes.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

load(fullfile(paths.data, 'wing_geom.mat'), 'cfg');

fprintf('=== AGLAS Stage 2: structural model ===\n\n');

fem   = beam_fem(cfg);
modes = modal_analysis(fem, min(12, fem.n_dof));

fprintf('Model: %d elements, %d DOF per node, %d free DOF\n', ...
        fem.n_el, fem.dpn, fem.n_dof);

% --- Sanity checks that the original code never performed --------------
tol = 1e-9;
assert(norm(fem.K - fem.K.', 'fro') < tol*norm(fem.K, 'fro'), 'K is not symmetric');
assert(norm(fem.M - fem.M.', 'fro') < tol*norm(fem.M, 'fro'), 'M is not symmetric');
assert(all(eig(fem.K) > 0), 'K is not positive definite');
assert(all(eig(fem.M) > 0), 'M is not positive definite');

mass_fem = sum(sum(fem.M_full(1:fem.dpn:end, 1:fem.dpn:end)));
mass_exact = cfg.section.m_total * cfg.geom.L;
fprintf('Mass check: assembled %.6f kg vs exact %.6f kg (error %.2e)\n', ...
        mass_fem, mass_exact, abs(mass_fem-mass_exact)/mass_exact);

orth_M = norm(modes.Phi.'*fem.M*modes.Phi - eye(modes.n_modes), 'fro');
fprintf('Mass-orthonormality  ||Phi''*M*Phi - I|| = %.3e\n', orth_M);

% --- Validation against the analytical cantilever ----------------------
beta_n = [1.875104068, 4.694091133, 7.854757438, 10.99554073, 14.13716839];
f_exact = beta_n.^2 / (2*pi) * ...
          sqrt(cfg.section.EI / (cfg.section.m_total * cfg.geom.L^4));
bend_idx = find(strcmp(modes.type, 'bending'));
n_cmp = min(numel(bend_idx), numel(f_exact));

fprintf('\nBending frequencies vs closed-form Euler-Bernoulli cantilever\n');
fprintf('  %-6s %12s %12s %10s\n', 'mode', 'FEM [Hz]', 'exact [Hz]', 'error [%]');
max_err = 0;
for i = 1:n_cmp
    f_fem = modes.freq_hz(bend_idx(i));
    err = 100*abs(f_fem - f_exact(i))/f_exact(i);
    max_err = max(max_err, err);
    fprintf('  %-6d %12.4f %12.4f %10.4f\n', i, f_fem, f_exact(i), err);
end
if max_err > 1.0
    warning('AGLAS:femAccuracy', ...
        'Bending frequency error %.2f %% exceeds 1 %%. Refine cfg.fem.n_elements.', max_err);
end

fprintf('\nAll modes\n');
for i = 1:min(8, modes.n_modes)
    fprintf('  mode %2d : %9.3f Hz   %s\n', i, modes.freq_hz(i), modes.type{i});
end

% --- Plot mode shapes ---------------------------------------------------
n_plot = min(4, modes.n_modes);
fig = figure('Name', 'Wing mode shapes', 'NumberTitle', 'off', ...
             'Position', [100 100 760 820]);
for i = 1:n_plot
    subplot(n_plot, 1, i);
    plot(fem.x_nodes, modes.w_shape(:,i), '-o', 'LineWidth', 1.6, 'MarkerSize', 4);
    if fem.torsion
        hold on;
        plot(fem.x_nodes, modes.ph_shape(:,i), '--s', 'LineWidth', 1.3, 'MarkerSize', 4);
        legend('bending w', 'twist \phi', 'Location', 'northwest');
    end
    grid on;
    ylabel('mode shape');
    title(sprintf('Mode %d  -  %.3f Hz  (%s)', i, modes.freq_hz(i), modes.type{i}));
    if i == n_plot, xlabel('spanwise station [m]'); end
end
save_figure(fig, fullfile(paths.results, 'mode_shapes.png'), 150);

out_file = fullfile(paths.data, 'structural_model.mat');
save(out_file, 'cfg', 'fem', 'modes');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'mode_shapes.png'));
