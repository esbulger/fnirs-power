function out = export_channel_info_NIRx(probe_file, montage_name, varargin)
%FNIRSPOWER.HELPERS.EXPORT_CHANNEL_INFO_NIRX
% Export NIRx channel names and positions for the Python power pipeline.
%
% Usage:
%   out = fnirspower.helpers.export_channel_info_NIRx( ...
%       probe_file, montage_name);
%
% Required inputs:
%   probe_file :
%       MAT file containing a variable named probeInfo.
%
%   montage_name :
%       Montage name used in the output filenames.
%
% Name-value options:
%   detector_cutoff :
%       Maximum detector index to retain. Use Inf to retain all channels.
%       Default: Inf.
%
%   output_dir :
%       Output directory. Default: <project>/python/data.
%
% Output:
%   out :
%       Struct containing channel_names, channel_positions, and the saved
%       file paths.
%
% Files written:
%   <montage_name>_channel_names.mat
%   <montage_name>_channel_positions.mat

P = fnirspower.paths();
default_output_dir = fullfile(P.root, 'python', 'data');

ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'probe_file', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addRequired(ip, 'montage_name', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

addParameter(ip, 'detector_cutoff', Inf, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(ip, 'output_dir', default_output_dir, ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));

parse(ip, probe_file, montage_name, varargin{:});
opts = ip.Results;

probe_file = char(opts.probe_file);
montage_name = matlab.lang.makeValidName(char(opts.montage_name));
output_dir = char(opts.output_dir);

if exist(probe_file, 'file') ~= 2
    error('Probe information file not found: %s', probe_file);
end

S = load(probe_file, 'probeInfo');

if ~isfield(S, 'probeInfo') || ~isfield(S.probeInfo, 'probes')
    error('Probe file must contain probeInfo.probes.');
end

probes = S.probeInfo.probes;

channel_mask = probes.index_c(:, 2) <= opts.detector_cutoff;
channel_link = probes.index_c(channel_mask, :);

if isempty(channel_link)
    error('No channels satisfied the requested detector cutoff.');
end

channel_name_cells = cell(size(channel_link, 1), 1);

for i = 1:size(channel_link, 1)
    channel_name_cells{i} = sprintf('%s-%s', ...
        probes.labels_s{channel_link(i, 1)}, ...
        probes.labels_d{channel_link(i, 2)});
end

channel_names = char(channel_name_cells);
channel_positions = probes.coords_c2(channel_mask, :);

if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

names_file = fullfile( ...
    output_dir, ...
    [montage_name '_channel_names.mat']);

positions_file = fullfile( ...
    output_dir, ...
    [montage_name '_channel_positions.mat']);

save(names_file, 'channel_names');
save(positions_file, 'channel_positions');

out = struct();
out.channel_names = channel_names;
out.channel_positions = channel_positions;
out.names_file = names_file;
out.positions_file = positions_file;
end
