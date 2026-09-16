% =========================================================================
% AGLAS SUPPLEMENTARY STAGE 2 : MESH CONVERGENCE STUDY
% =========================================================================
% Establishes that the finite element discretisation is converged, and
% separates numerical error from physical modelling effects.
%
% Three checks:
%   A  Pure-bending frequencies against the closed-form Euler-Bernoulli
%      cantilever. The torsion DOF is switched off here so the FE model and
%      the analytical reference describe the same physics.
%   B  The frequency shift produced by inertial bending-torsion coupling,
%      which is a physical effect of the CG being offset from the elastic
%      axis, not a discretisation error. Reported separately so it is not
%      mistaken for one.
%   C  Static response: tip deflection against q*L^4/(8*EI) and root bending
%      moment against q*L^2/2.
%
% This is the evidence for the element count used elsewhere. The original
% code chose 10 elements with no convergence study at all.
%
% Outputs : data/mesh_convergence.mat, results/mesh_convergence.png
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

cfg = aglas_config('baseline');

fprintf('=== AGLAS supplementary stage 2: mesh convergence ===\n\n');

n_el_list = [2 4 5 8 10 16 20 32 40 64];
n_cases   = numel(n_el_list);

beta_n  = [1.875104068, 4.694091133, 7.854757438];
f_exact = beta_n.^2 / (2*pi) * ...
          sqrt(cfg.section.EI / (cfg.section.m_total * cfg.geom.L^4));

q_ref   = 1000;                                   % reference line load [N/m]
w_exact = q_ref * cfg.geom.L^4 / (8 * cfg.section.EI);
M_exact = q_ref * cfg.geom.L^2 / 2;

err_f  = nan(n_cases, 3);     % pure bending model
err_w  = nan(n_cases, 1);
err_M  = nan(n_cases, 1);
f_bend = nan(n_cases, 3);

% ---------------------------------------------------- A: pure bending --
fprintf('A. Pure-bending model vs closed-form cantilever (torsion disabled)\n');
fprintf('%9s %11s %11s %11s %13s %12s\n', ...
        'elements', 'e_f1 [%]', 'e_f2 [%]', 'e_f3 [%]', 'e_tip [%]', 'e_Mroot [%]');

for i = 1:n_cases
    c = cfg;
    c.fem.n_elements = n_el_list(i);
    c.fem.include_torsion = false;

    fem = beam_fem(c);
    md  = modal_analysis(fem, min(12, fem.n_dof));

    nb = min(3, md.n_modes);
    f_bend(i,1:nb) = md.freq_hz(1:nb).';
    err_f(i,1:nb)  = 100*abs(f_bend(i,1:nb) - f_exact(1:nb)) ./ f_exact(1:nb);

    f_vec = consistent_line_load(fem, q_ref, 0);
    u     = fem.K \ f_vec;
    err_w(i) = 100*abs(u(fem.idx_w(end)) - w_exact)/w_exact;

    rec = recover_bending(c, fem, u);
    err_M(i) = 100*abs(rec.moment(1) - M_exact)/M_exact;

    fprintf('%9d %11.3e %11.3e %11.3e %13.2e %12.4f\n', ...
            n_el_list(i), err_f(i,1), err_f(i,2), err_f(i,3), err_w(i), err_M(i));
end

fit_range = n_el_list >= 4 & n_el_list <= 40;
p_f1 = polyfit(log(n_el_list(fit_range)), log(err_f(fit_range,1).'), 1);
p_M  = polyfit(log(n_el_list(fit_range)), log(err_M(fit_range).'), 1);

fprintf('\n  Observed order of convergence\n');
fprintf('    fundamental frequency : %.2f\n', -p_f1(1));
fprintf('    root bending moment   : %.2f   (theory: 2)\n', -p_M(1));
fprintf('  Tip deflection is exact to machine precision at every mesh density,\n');
fprintf('  because the cubic Hermite element reproduces the exact displacement\n');
fprintf('  field of a uniformly loaded prismatic beam.\n');

% --------------------------------------- B: bending-torsion coupling --
fprintf('\nB. Effect of inertial bending-torsion coupling (physical, not numerical)\n');

c_fine = cfg; c_fine.fem.n_elements = 64;
f_coupled = modal_analysis(beam_fem(c_fine), 6).freq_hz;

c_bal = c_fine;
c_bal.geom.cg_frac = c_bal.geom.ea_frac;       % CG moved onto the elastic axis
c_bal.aero.x_alpha = 0;
c_bal.section = section_properties(c_bal);
f_balanced = modal_analysis(beam_fem(c_bal), 6).freq_hz;

fprintf('  CG offset from elastic axis : %.1f %% chord (%.4f m)\n', ...
        100*(cfg.geom.cg_frac - cfg.geom.ea_frac), cfg.section.d_ea_cg);
fprintf('  first bending, CG offset    : %.6f Hz\n', f_coupled(1));
fprintf('  first bending, mass balanced: %.6f Hz\n', f_balanced(1));
fprintf('  shift                       : %.4f %%\n', ...
        100*(f_coupled(1) - f_balanced(1))/f_balanced(1));
fprintf('  This shift is why the coupled model does not converge onto the\n');
fprintf('  pure-bending analytical value: the two describe different systems.\n');
fprintf('  It is also the coupling that makes flutter possible in stage 4.\n');

% --------------------------------------------------------- conclusion --
k = find(n_el_list == cfg.fem.n_elements, 1);
if ~isempty(k)
    fprintf('\nAt the configured %d elements: f1 error %.2e %%, root moment error %.4f %%\n', ...
            cfg.fem.n_elements, err_f(k,1), err_M(k));
end

% ------------------------------------------------------------- plots --
fig = figure('Name', 'Mesh convergence', 'NumberTitle', 'off', ...
             'Position', [100 100 950 420]);

subplot(1,2,1);
loglog(n_el_list, err_f(:,1), '-o', 'LineWidth', 1.6); hold on;
loglog(n_el_list, err_f(:,2), '-s', 'LineWidth', 1.4);
loglog(n_el_list, err_f(:,3), '-^', 'LineWidth', 1.4);
grid on; xlabel('number of elements'); ylabel('frequency error [%]');
title('Natural frequency convergence (pure bending)');
legend('mode 1', 'mode 2', 'mode 3', 'Location', 'southwest');

subplot(1,2,2);
loglog(n_el_list, err_M, '-o', 'LineWidth', 1.6); hold on;
ref = err_M(2) * (n_el_list/n_el_list(2)).^(-2);
loglog(n_el_list, ref, 'k--', 'LineWidth', 1.1);
grid on; xlabel('number of elements'); ylabel('root moment error [%]');
title('Recovered root bending moment convergence');
legend('measured', 'second order reference', 'Location', 'southwest');

save_figure(fig, fullfile(paths.results, 'mesh_convergence.png'), 150);

out_file = fullfile(paths.data, 'mesh_convergence.mat');
save(out_file, 'cfg', 'n_el_list', 'err_f', 'err_w', 'err_M', 'f_bend', ...
     'f_exact', 'f_coupled', 'f_balanced');
fprintf('\nSaved %s\n', out_file);
fprintf('Saved %s\n', fullfile(paths.results, 'mesh_convergence.png'));
