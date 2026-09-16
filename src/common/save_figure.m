function save_figure(fig, filename, dpi)
%SAVE_FIGURE  Export a figure to PNG with a light background.
%
%   SAVE_FIGURE(fig, filename)
%   SAVE_FIGURE(fig, filename, dpi)
%
%   Wraps print so every stage exports the same way, and so the result does not
%   depend on whether the person running the pipeline has the MATLAB desktop in
%   light or dark mode. See FORCE_LIGHT_FIGURE.

    if nargin < 3 || isempty(dpi), dpi = 150; end

    force_light_figure(fig);
    print(fig, filename, '-dpng', sprintf('-r%d', dpi));
end
