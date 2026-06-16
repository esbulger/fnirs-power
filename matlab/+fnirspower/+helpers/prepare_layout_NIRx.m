function layout = prepare_layout_NIRx(probe_file, montage_name, varargin)
%FNIRSPOWER.HELPERS.PREPARE_LAYOUT_NIRX
% Prepare and optionally save a 2D FieldTrip layout from NIRx probe data.
%
% Usage:
%   layout_dense = fnirspower.helpers.prepare_layout_NIRx( ...
%       probe_file, 'dense');
%
% Required inputs:
%   probe_file    Path to a MAT file containing a variable named probeInfo.
%                 Psrent folder can be in NIRx/configurations directory.
%
%   montage_name  Name used to construct the saved filename and variable:
%                 layout_<montage_name>.
%
% Name-value options:
%   output_file              MAT-file path used to save the layout.
%                            Default:
%                            fnirspower/workspace/layouts/layout_<montage_name>.mat
%                            Use [] to disable saving.
%
%   detector_cutoff          Retain detector indices below this value.
%                            Default: 999.
%
%   position_scale           Multiplicative [x, y] position scaling.
%                            Default: [1.1, 1.1].
%
%   position_shift           Additive [x, y] position shift.
%                            Default: [0, 0].
%
%   layout_scaling           Multiplicative scaling of channel marker
%                            width and height. Default: 1.
%
%   build_hemisphere_masks   Replace the FieldTrip mask with separate
%                            left/right convex-hull masks. Default: false.
%
%   mask_offset              Outward mask expansion in layout units.
%                            Default: 0.025.
%
%   show_demo                Plot simulated channel data. Default: false.
%
% Output:
%   layout   FieldTrip layout struct.
%
% Notes:
%   MATLAB function outputs cannot be dynamically named in the caller's
%   workspace. The caller chooses that name on the left-hand side. When
%   saved, however, both the MAT filename and the variable stored inside
%   the file are named layout_<montage_name>.

if nargin < 2
    error('prepare_layout_NIRx requires probe_file and montage_name.');
end

if ~(ischar(probe_file) || ...
        (isstring(probe_file) && isscalar(probe_file)))
    error('probe_file must be a character vector or string scalar.');
end

if ~(ischar(montage_name) || ...
        (isstring(montage_name) && isscalar(montage_name)))
    error('montage_name must be a character vector or string scalar.');
end

probe_file = char(probe_file);
montage_name = strtrim(char(montage_name));

if isempty(montage_name)
    error('montage_name cannot be empty.');
end

% Construct a valid MATLAB variable name and matching default filename.
layout_name = matlab.lang.makeValidName(['layout_' montage_name]);

P = fnirspower.paths();
default_output_file = fullfile(P.layouts, [layout_name '.mat']);

ip = inputParser;
ip.FunctionName = mfilename;

addParameter(ip, 'output_file', default_output_file, ...
    @(x) isempty(x) || ischar(x) || ...
    (isstring(x) && isscalar(x)));

addParameter(ip, 'detector_cutoff', 999, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

addParameter(ip, 'position_scale', [1.1, 1.1], ...
    @(x) isnumeric(x) && numel(x) == 2 && ...
    all(isfinite(x)) && all(x > 0));

addParameter(ip, 'position_shift', [0, 0], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)));

addParameter(ip, 'layout_scaling', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

addParameter(ip, 'build_hemisphere_masks', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

addParameter(ip, 'mask_offset', 0.025, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);

addParameter(ip, 'show_demo', false, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

parse(ip, varargin{:});
opts = ip.Results;

opts.position_scale = reshape(opts.position_scale, 1, 2);
opts.position_shift = reshape(opts.position_shift, 1, 2);

if exist('ft_prepare_layout', 'file') ~= 2
    error('FieldTrip must be available on the MATLAB path.');
end

if exist(probe_file, 'file') ~= 2
    error('Probe information file not found: %s', probe_file);
end

loaded = load(probe_file, 'probeInfo');

if ~isfield(loaded, 'probeInfo')
    error('The probe MAT file must contain a variable named probeInfo.');
end

probeInfo = loaded.probeInfo;

if ~isfield(probeInfo, 'probes') || ~isstruct(probeInfo.probes)
    error('probeInfo must contain a probes struct.');
end

probes = probeInfo.probes;
required_fields = {'index_c', 'coords_c3', 'labels_s', 'labels_d'};

for k = 1:numel(required_fields)
    if ~isfield(probes, required_fields{k})
        error('probeInfo.probes is missing the field "%s".', ...
            required_fields{k});
    end
end

% Retain long channels using the montage-specific detector rule.
channel_mask = probes.index_c(:, 2) < opts.detector_cutoff;
long_chan_link = probes.index_c(channel_mask, :);
n_channels_layout = size(long_chan_link, 1);

if n_channels_layout == 0
    error('No channels satisfied detector index < %g.', ...
        opts.detector_cutoff);
end

if size(probes.coords_c3, 1) ~= size(probes.index_c, 1)
    error(['probes.coords_c3 must contain one coordinate row for each ', ...
        'row of probes.index_c.']);
end

labels = cell(n_channels_layout, 1);

for channel_idx = 1:n_channels_layout
    source_idx = long_chan_link(channel_idx, 1);
    detector_idx = long_chan_link(channel_idx, 2);

    labels{channel_idx} = sprintf('%s-%s', ...
        local_get_label(probes.labels_s, source_idx), ...
        local_get_label(probes.labels_d, detector_idx));
end

opto = struct();
opto.label = labels;
opto.chanpos = probes.coords_c3(channel_mask, :);
opto.elecpos = probes.coords_c3(channel_mask, :);
opto.unit = 'mm';

cfg = [];
cfg.elec = opto;
layout = ft_prepare_layout(cfg);

% Remove FieldTrip's optional COMNT/SCALE entries.
layout.pos = layout.pos(1:n_channels_layout, :);
layout.label = layout.label(1:n_channels_layout);
layout.width = layout.width(1:n_channels_layout);
layout.height = layout.height(1:n_channels_layout);

% Apply configurable position and marker adjustments.
layout.pos(:, 1:2) = ...
    layout.pos(:, 1:2) .* opts.position_scale;

layout.pos(:, 1:2) = ...
    layout.pos(:, 1:2) + opts.position_shift;

layout.width = layout.width * opts.layout_scaling;
layout.height = layout.height * opts.layout_scaling;

% Retain FieldTrip's default mask unless custom hemisphere masks are
% explicitly requested.
if logical(opts.build_hemisphere_masks)
    channel_positions = layout.pos(:, 1:2);
    positions_left = channel_positions(channel_positions(:, 1) < 0, :);
    positions_right = channel_positions(channel_positions(:, 1) >= 0, :);

    if size(positions_left, 1) < 3 || size(positions_right, 1) < 3
        error(['At least three channels per side are required to build ', ...
            'hemisphere masks.']);
    end

    layout.mask = { ...
        local_offset_hull(positions_left, opts.mask_offset); ...
        local_offset_hull(positions_right, opts.mask_offset)};
end

% Save both the MAT filename and contained variable using layout_name.
if ~isempty(opts.output_file)
    output_file = char(opts.output_file);
    output_dir = fileparts(output_file);

    if ~isempty(output_dir) && exist(output_dir, 'dir') ~= 7
        mkdir(output_dir);
    end

    save_data = struct();
    save_data.(layout_name) = layout;
    save(output_file, '-struct', 'save_data');
end

if logical(opts.show_demo)
    local_plot_demo(layout);
end
end


function label = local_get_label(label_array, index)
% Return one source or detector label as a character vector.

if index < 1 || index > numel(label_array)
    error('Probe label index %d is out of range.', index);
end

if iscell(label_array)
    label = label_array{index};
elseif isstring(label_array)
    label = label_array(index);
elseif ischar(label_array)
    label = label_array(index, :);
else
    error('Probe labels must be a cell, string, or character array.');
end

label = strtrim(char(label));
end


function polygon = local_offset_hull(points, offset)
% Create a closed convex-hull polygon with a small radial expansion.

hull_idx = convhull(points(:, 1), points(:, 2));
polygon = points(hull_idx, :);

if offset > 0
    center = mean(polygon(1:end-1, :), 1);
    direction = polygon - center;
    distance = sqrt(sum(direction.^2, 2));
    polygon = polygon + offset .* direction ./ max(distance, eps);
end
end


function local_plot_demo(layout)
% Plot simulated channel data using the prepared layout.

if exist('ft_topoplotER', 'file') ~= 2
    error('ft_topoplotER must be available to display the demo plot.');
end

cfg = [];
cfg.layout = layout;
cfg.box = 'no';
cfg.comment = 'no';
cfg.parameter = 'avg';
cfg.interpolation = 'v4';

timelock = [];
timelock.avg = randn(numel(layout.label), 1);
timelock.time = 1;
timelock.label = layout.label;
timelock.dimord = 'chan_time';

figure;
ft_topoplotER(cfg, timelock);
colorbar;
caxis([-5, 5]);
end
