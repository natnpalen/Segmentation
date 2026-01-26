function scroll_viewer(volume_data, ds, title_prefix, cmap, clabel, varargin)
% SCROLL_VIEWER  Triplanar scroll viewer (migrated from run_density_analysis helpers).

create_interactive_viewer(volume_data, 3, ds, [title_prefix, ' - Axial'],    cmap, clabel, varargin{:});
create_interactive_viewer(volume_data, 2, ds, [title_prefix, ' - Sagittal'], cmap, clabel, varargin{:});
create_interactive_viewer(volume_data, 1, ds, [title_prefix, ' - Coronal'],  cmap, clabel, varargin{:});

end

function create_interactive_viewer(volume, view_dim, ds, fig_title, cmap, clabel, varargin)
p = inputParser;
addParameter(p,'clim',[],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'alpha',[],@isnumeric);
parse(p,varargin{:});
opts = p.Results;

sz = size(volume);
num_slices  = sz(view_dim);
curSlice    = max(1, min(num_slices, round(num_slices/2)));

fig = figure('Name', fig_title, 'Color', 'w', 'Toolbar', 'figure');
ax  = axes('Parent', fig);

S0 = get_slice(volume, curSlice, view_dim);
hImg = imagesc(ax, S0);
if ~isempty(opts.alpha)
    A0 = get_slice(opts.alpha, curSlice, view_dim);
    set(hImg, 'AlphaData', A0);
end
axis(ax,'tight'); axis(ax,'off');
colormap(ax, cmap);
c = colorbar(ax); c.Label.String = clabel;
if ~isempty(opts.clim)
    clim(ax, opts.clim);
else
    num_data_colors = size(cmap,1) - 1;
    if any(volume(:) > num_data_colors) || any(volume(:) < 0)
        % autoscale
    else
        clim(ax, [0 num_data_colors]);
    end
end
drawnow;

slider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units','normalized', 'Position', [0.25 0.02 0.50 0.04], ...
    'Min', 1, 'Max', num_slices, 'Value', curSlice, ...
    'SliderStep', [1/max(1,num_slices-1), 10/max(1,num_slices-1)], ...
    'Interruptible','off','BusyAction','cancel', ...
    'Callback', @(src,~) setSlice(round(src.Value)));

set(fig, 'WindowScrollWheelFcn', @(~,evt) onScroll(evt));
set(fig, 'KeyPressFcn',        @(~,evt) onKey(evt));
update_title(ax, curSlice, view_dim);

    function setSlice(n)
        n = max(1, min(num_slices, n));
        if n == curSlice, return; end
        curSlice = n;
        hImg.CData = get_slice(volume, curSlice, view_dim);
        if ~isempty(opts.alpha)
            hImg.AlphaData = get_slice(opts.alpha, curSlice, view_dim);
        end
        slider.Value = curSlice;
        update_title(ax, curSlice, view_dim);
        drawnow limitrate;
    end
    function onScroll(evt)
        step = 1;
        if ismember('shift', get(fig,'CurrentModifier')), step = 5; end
        setSlice(curSlice + sign(evt.VerticalScrollCount)*(-step));
    end
    function onKey(evt)
        switch evt.Key
            case {'rightarrow','uparrow'},   setSlice(curSlice+1);
            case {'leftarrow','downarrow'},  setSlice(curSlice-1);
            case 'pageup',                   setSlice(curSlice+5);
            case 'pagedown',                 setSlice(curSlice-5);
            case 'home',                     setSlice(1);
            case 'end',                      setSlice(num_slices);
        end
    end
end

function slice_data = get_slice(volume, slice_num, view_dim)
slice_num = round(slice_num);
switch view_dim
    case 1, slice_data = squeeze(volume(slice_num, :, :))';
    case 2, slice_data = squeeze(volume(:, slice_num, :))';
    case 3, slice_data = squeeze(volume(:, :, slice_num));
end
end

function update_title(ax, slice_num, view_dim)
slice_num = round(slice_num);
switch view_dim
    case 1, title(ax, sprintf('Coronal (Row %d)', slice_num));
    case 2, title(ax, sprintf('Sagittal (Column %d)', slice_num));
    case 3, title(ax, sprintf('Axial (Slice %d)', slice_num));
end
end
