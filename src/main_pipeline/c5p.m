% =========================================================================
% AGLAS STAGE 5 (main pipeline) : DEFLECTION ANIMATION
% =========================================================================
% Animates the reconstructed deflected shape over the gust event and writes
% an animated GIF to animation/.
%
% Inputs  : data/modal_response.mat
% Outputs : animation/wing_gust_animation.gif
% =========================================================================

clear; close all; clc;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'common'));
paths = aglas_paths();

SAVE_GIF    = true;
MAX_FRAMES  = 80;      % cap the frame count so the GIF stays a sane size
FRAME_DELAY = 0.06;    % seconds between GIF frames

load(fullfile(paths.data, 'modal_response.mat'), ...
     'cfg', 't_sol', 'w_field', 'tw', 'wg', 't_grid');

fprintf('=== AGLAS Stage 5: animation ===\n\n');

x        = linspace(0, cfg.geom.L, size(w_field, 2));
frame_step = max(1, ceil(numel(t_sol)/MAX_FRAMES));
frames     = 1:frame_step:numel(t_sol);
gif_file = fullfile(paths.animation, 'wing_gust_animation.gif');

y_lim = max(abs(w_field(:))) * 1.25;
if y_lim == 0 || ~isfinite(y_lim), y_lim = 1; end

fig = figure('Name', 'Wing gust response', 'NumberTitle', 'off', ...
             'Position', [100 100 720 430], 'Color', 'w');

fprintf('Rendering %d frames ...\n', numel(frames));
first = true;
for idx = frames
    clf(fig);

    subplot(2,1,1);
    plot(x, w_field(idx,:), '-', 'LineWidth', 2.4); hold on;
    plot(x, zeros(size(x)), 'k--', 'LineWidth', 0.8);
    grid on;
    xlim([0 cfg.geom.L]); ylim([-y_lim y_lim]);
    xlabel('spanwise station [m]'); ylabel('deflection [m]');
    title(sprintf('t = %5.2f s     tip deflection = %+6.3f m', ...
                  t_sol(idx), w_field(idx,end)));

    subplot(2,1,2);
    plot(t_grid, wg, 'LineWidth', 1.4); hold on;
    plot(t_sol(idx), interp1(t_grid, wg, t_sol(idx)), 'o', ...
         'MarkerSize', 8, 'LineWidth', 1.6);
    grid on;
    xlim([t_grid(1) t_grid(end)]);
    xlabel('time [s]'); ylabel('gust velocity [m/s]');
    title('Gust input');

    % Reapply after every clf: the cleared figure creates fresh axes, which
    % pick up the MATLAB desktop theme again. Without this the GIF comes out
    % black-backgrounded for anyone running in dark mode.
    force_light_figure(fig);

    drawnow;

    if SAVE_GIF
        frame = getframe(fig);
        [ind, cmap] = rgb_to_indexed(frame2im(frame), 128);
        if first
            imwrite(ind, cmap, gif_file, 'gif', ...
                    'LoopCount', Inf, 'DelayTime', FRAME_DELAY);
            first = false;
        else
            imwrite(ind, cmap, gif_file, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', FRAME_DELAY);
        end
    else
        pause(FRAME_DELAY);
    end
end

if SAVE_GIF
    fprintf('Saved %s\n', gif_file);
end

