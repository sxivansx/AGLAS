function force_light_figure(fig)
%FORCE_LIGHT_FIGURE  Give a figure a light background whatever the desktop theme.
%
%   FORCE_LIGHT_FIGURE(fig)
%
%   From R2025a, MATLAB has a dark desktop theme and figures inherit it. print
%   then exports exactly what is on screen, so anyone running with dark mode on
%   gets black-backgrounded PNGs. That is fine on screen and wrong for a report,
%   a LaTeX document, or a README rendered on a white page.
%
%   This forces the figure chrome to a light scheme before export: figure and
%   axes backgrounds white, axis lines, ticks, titles and labels black, legends
%   and colorbars to match. Plotted data colours are left alone, as are text
%   objects the caller coloured deliberately, such as the red flutter marker.
%
%   Every property is set behind an isprop guard, so the function degrades to a
%   no-op on releases and interpreters that lack a given property.

    if nargin < 1 || isempty(fig), fig = gcf; end

    black     = [0 0 0];
    grid_grey = [0.15 0.15 0.15];

    % R2025a and later expose an explicit theme switch.
    try %#ok<TRYNC>
        theme(fig, 'light');
    end

    set(fig, 'Color', 'w');
    try %#ok<TRYNC>
        set(fig, 'InvertHardcopy', 'on');
    end

    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        a = ax(k);
        set_if(a, 'Color',          'w');
        set_if(a, 'XColor',         black);
        set_if(a, 'YColor',         black);
        set_if(a, 'ZColor',         black);
        set_if(a, 'GridColor',      grid_grey);
        set_if(a, 'MinorGridColor', grid_grey);
        for p = {'Title', 'XLabel', 'YLabel', 'ZLabel'}
            if isprop(a, p{1})
                h = get(a, p{1});
                if ~isempty(h) && all(ishandle(h))
                    set_if(h, 'Color', black);
                end
            end
        end
    end

    lg = find_of_type(fig, 'legend');
    for k = 1:numel(lg)
        set_if(lg(k), 'Color',     'w');
        set_if(lg(k), 'TextColor', black);
        set_if(lg(k), 'EdgeColor', grid_grey);
    end

    % A colorbar's Color property is its tick and label colour, not a
    % background, so it takes black rather than white.
    cb = find_of_type(fig, 'colorbar');
    for k = 1:numel(cb)
        set_if(cb(k), 'Color', black);
    end
end

% ------------------------------------------------------------------------
function h = find_of_type(fig, type_name)
%FIND_OF_TYPE  findall by type, tolerant of interpreters that lack the type.
    h = [];
    try %#ok<TRYNC>
        h = findall(fig, 'Type', type_name);
    end
end

% ------------------------------------------------------------------------
function set_if(h, name, value)
%SET_IF  Set a property only when the handle actually has it.
    try %#ok<TRYNC>
        if isprop(h, name)
            set(h, name, value);
        end
    end
end
