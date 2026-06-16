function fig_handles = plot_channel_level_hb( ...
    layout, data2use, clim, date_value, task_type, montage_type, ...
    hbr, interp_limit, output_dir)
%FNIRSPOWER.HELPERS.PLOT_CHANNEL_LEVEL_HB
% Plot channel-level HbO and optional HbR topographies using FieldTrip.
%
% Usage:
%   fig_handles = fnirspower.helpers.plot_channel_level_hb( ...
%       layout, data2use, clim, date_value, task_type, montage_type, hbr)
%
%   fig_handles = fnirspower.helpers.plot_channel_level_hb( ...
%       layout, data2use, clim, date_value, task_type, montage_type, hbr, ...
%       interp_limit, output_dir)
%
% Inputs:
%   layout         FieldTrip layout struct containing channel labels and
%                  positions.
%
%   data2use       Channel-level values [1 x nChannels] or
%                  [2 x nChannels]. Row 1 contains HbO; row 2 contains HbR.
%
%   clim           Positive scalar defining symmetric color limits:
%                  [-clim, clim].
%
%   date_value     Date text used in output filenames. If empty, the
%                  current date is used.
%
%   task_type      Task label used in plot titles and output filenames.
%
%   montage_type   Montage label used in output filenames.
%
%   hbr            If true, also plot data2use row 2 as HbR.
%
%   interp_limit   Optional FieldTrip interplimits value. Default:
%                  'sensors'.
%
%   output_dir     Optional output directory. Default:
%                  <project>/figures
%                  Pass [] to display figures without saving them.
%
% Output:
%   fig_handles    Struct containing figure handles:
%                    .hbo
%                    .hbr
%
% Dependencies:
%   FieldTrip functions ft_topoplotER and ft_colormap.

if nargin < 8 || isempty(interp_limit)
    interp_limit = 'sensors';
end

if nargin < 9
    P = fnirspower.paths();
    output_dir = P.figures;
end

if exist('ft_topoplotER', 'file') ~= 2
    error('FieldTrip function ft_topoplotER was not found on the MATLAB path.');
end

if ~isstruct(layout) || ~isfield(layout, 'label')
    error('layout must be a FieldTrip layout struct containing layout.label.');
end

if ~isnumeric(data2use) || isempty(data2use) || ndims(data2use) > 2
    error('data2use must be a nonempty numeric matrix.');
end

n_channels = numel(layout.label);

if isvector(data2use) && numel(data2use) == n_channels
    data2use = reshape(data2use, 1, n_channels);
end

if size(data2use, 2) ~= n_channels
    error(['data2use must contain one column per layout channel. ', ...
        'Expected %d columns but received %d.'], ...
        n_channels, size(data2use, 2));
end

if ~(isnumeric(clim) && isscalar(clim) && isfinite(clim) && clim > 0)
    error('clim must be a finite positive scalar.');
end

if ~(islogical(hbr) || (isnumeric(hbr) && isscalar(hbr)))
    error('hbr must be a logical or numeric scalar.');
end

plot_hbr = logical(hbr);

if plot_hbr && size(data2use, 1) < 2
    error('HbR plotting requires a second row in data2use.');
end

if ~(ischar(interp_limit) || ...
        (isstring(interp_limit) && isscalar(interp_limit)))
    error('interp_limit must be a character vector or string scalar.');
end

if ~(isempty(output_dir) || ischar(output_dir) || ...
        (isstring(output_dir) && isscalar(output_dir)))
    error('output_dir must be empty, a character vector, or a string scalar.');
end

if isempty(date_value)
    date_text = char(string(datetime('now', 'Format', 'yyyy-MM-dd')));
else
    date_text = char(string(date_value));
end

task_text = char(string(task_type));
montage_text = char(string(montage_type));

date_file = local_filename_component(date_text);
task_file = local_filename_component(task_text);
montage_file = local_filename_component(montage_text);

if ~isempty(output_dir)
    output_dir = char(output_dir);

    if exist(output_dir, 'dir') ~= 7
        mkdir(output_dir);
    end
end

% Shared FieldTrip plotting configuration.
cfg = [];
cfg.layout = layout;
cfg.box = 'no';
cfg.comment = 'no';
cfg.parameter = 'avg';
cfg.interpolation = 'v4';
cfg.interplimits = char(interp_limit);
cfg.zlim = [-clim, clim];
cfg.maskparameter = 'nan';

fig_handles = struct();
fig_handles.hbo = gobjects(1);
fig_handles.hbr = gobjects(0);

% Plot HbO.
timelock = local_make_timelock(data2use(1, :), layout.label);

fig_handles.hbo = figure;
ft_topoplotER(cfg, timelock);

colorbar_handle = colorbar;
caxis([-clim, clim]);
ft_colormap('redblue');

title(['HbO_' task_text], 'Interpreter', 'none');
colorbar_handle.Label.String = '\DeltaHbO (\muM)';

if ~isempty(output_dir)
    output_file = fullfile( ...
        output_dir, ...
        sprintf('%s_HbO_Label_%s_%s.png', ...
            date_file, task_file, montage_file));

    exportgraphics(gca, output_file);
end

% Plot HbR when requested.
if plot_hbr
    timelock = local_make_timelock(data2use(2, :), layout.label);

    fig_handles.hbr = figure;
    ft_topoplotER(cfg, timelock);

    colorbar_handle = colorbar;
    caxis([-clim, clim]);
    ft_colormap('redblue');

    title(['HbR_' task_text], 'Interpreter', 'none');
    colorbar_handle.Label.String = '\DeltaHbR (\muM)';

    if ~isempty(output_dir)
        output_file = fullfile( ...
            output_dir, ...
            sprintf('%s_HbR_Label_%s_%s.png', ...
                date_file, task_file, montage_file));

        exportgraphics(gca, output_file);
    end
end
end


function timelock = local_make_timelock(channel_values, labels)
% Construct a single-time-point FieldTrip timelock structure.

timelock = [];
timelock.avg = channel_values(:);
timelock.time = 1;
timelock.label = labels;
timelock.dimord = 'chan_time';
end


function value = local_filename_component(value)
% Replace filename-unsafe characters while preserving readable labels.

value = strtrim(char(value));
value = regexprep(value, '[^\w.-]+', '_');

if isempty(value)
    value = 'unspecified';
end
end
