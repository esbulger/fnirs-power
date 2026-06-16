function [subj_prob, avg_prob, roi_info] = compute_subject_ROI( ...
    p, roi_def, varargin)
%FNIRSPOWER.MEASPRED.COMPUTE_SUBJECT_ROI
% Generate subject-specific spherical ROI maps on a brain mesh.
%
% For each subject and ROI:
%   1) Translate the nominal ROI center independently along x, y, and z.
%   2) Limit each translation to within two standard deviations for realistic modelling.
%   3) Snap the translated center to the nearest brain node.
%   4) Draw a subject-specific ROI radius.
%   5) Select all brain nodes within the resulting sphere.
%
% ROI definition fields
% ---------------------
% roi_def.center or roi_def.center_mm :
%     Single ROI center [1 x 3] or [3 x 1].
%
% roi_def.centers or roi_def.centers_mm :
%     Multiple ROI centers [nROI x 3].
%
% roi_def.radius_mm :
%     Mean spherical ROI radius in millimeters. May be scalar or contain
%     one value per ROI. Default: 2.
%
% roi_def.radius_std_mm :
%     Standard deviation of the subject-specific ROI radius in
%     millimeters. May be scalar or contain one value per ROI.
%     Default: 0.
%
% roi_def.sigma_mm :
%     Standard deviation of the subject-specific center translation along
%     each coordinate axis, in millimeters. May be scalar or contain one
%     value per ROI. Default: 5.
%
% Options
% -------
% n_subj :
%     Number of subject-specific ROI maps. Default: 16.
%
% seed :
%     Optional random seed. Default: [].
%
% Outputs
% -------
% subj_prob :
%     Binary subject-specific ROI maps
%     [nSubjects x nROI x nNodes].
%
% avg_prob :
%     Mean ROI inclusion probability across subjects
%     [nROI x nNodes].
%
% roi_info :
%     Struct containing the nominal centers, translated centers, nearest
%     mesh-node centers, subject-specific radii, and resolved options.

% -------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------
if nargin < 2
    error('compute_subject_ROI requires p and roi_def.');
end

if isempty(p) || ~isnumeric(p) || size(p, 2) ~= 3
    error('p must contain nNodes-by-3 numeric node coordinates.');
end

if ~isstruct(roi_def) || ~isscalar(roi_def)
    error('roi_def must be a scalar struct.');
end

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
opts = struct();
opts.n_subj = 16;
opts.seed = [];

name_value_inputs = varargin;

if ~isempty(name_value_inputs) && isstruct(name_value_inputs{1})
    supplied_opts = name_value_inputs{1};
    name_value_inputs = name_value_inputs(2:end);

    supplied_fields = fieldnames(supplied_opts);

    for field_idx = 1:numel(supplied_fields)
        field_name = supplied_fields{field_idx};

        if ~isfield(opts, field_name)
            error('Unknown option "%s".', field_name);
        end

        opts.(field_name) = supplied_opts.(field_name);
    end
end

if ~isempty(name_value_inputs)
    if mod(numel(name_value_inputs), 2) ~= 0
        error('Name-value inputs must be supplied in pairs.');
    end

    for input_idx = 1:2:numel(name_value_inputs)
        option_name = name_value_inputs{input_idx};
        option_value = name_value_inputs{input_idx + 1};

        if ~(ischar(option_name) || ...
                (isstring(option_name) && isscalar(option_name)))
            error('Option names must be character vectors or string scalars.');
        end

        option_name = char(option_name);

        if ~isfield(opts, option_name)
            error('Unknown option "%s".', option_name);
        end

        opts.(option_name) = option_value;
    end
end

if ~(isnumeric(opts.n_subj) && isscalar(opts.n_subj) && ...
        isfinite(opts.n_subj) && opts.n_subj >= 1 && ...
        opts.n_subj == floor(opts.n_subj))
    error('n_subj must be a positive integer.');
end

% Use a local random-number stream when a seed is supplied.
random_stream = [];

if ~isempty(opts.seed)
    random_stream = RandStream( ...
        'mt19937ar', ...
        'Seed', opts.seed);
end

% -------------------------------------------------------------------------
% Resolve ROI centers and parameters
% -------------------------------------------------------------------------
centers = local_resolve_centers(roi_def);

n_roi = size(centers, 1);
n_nodes = size(p, 1);
n_subjects = opts.n_subj;

radius_mm = local_get_roi_parameter( ...
    roi_def, ...
    'radius_mm', ...
    n_roi, ...
    2);

radius_std_mm = local_get_roi_parameter( ...
    roi_def, ...
    'radius_std_mm', ...
    n_roi, ...
    0);

sigma_mm = local_get_roi_parameter( ...
    roi_def, ...
    'sigma_mm', ...
    n_roi, ...
    5);

if any(~isfinite(radius_mm)) || any(radius_mm <= 0)
    error('roi_def.radius_mm values must be finite and greater than zero.');
end

if any(~isfinite(radius_std_mm)) || any(radius_std_mm < 0)
    error('roi_def.radius_std_mm values must be finite and nonnegative.');
end

if any(~isfinite(sigma_mm)) || any(sigma_mm < 0)
    error('roi_def.sigma_mm values must be finite and nonnegative.');
end

% -------------------------------------------------------------------------
% Generate subject-specific ROI maps
% -------------------------------------------------------------------------
subj_prob = zeros(n_subjects, n_roi, n_nodes);

center_translations = zeros(n_roi, 3, n_subjects);
translated_centers = NaN(n_roi, 3, n_subjects);
mesh_centers = NaN(n_roi, 3, n_subjects);
mesh_center_indices = NaN(n_roi, n_subjects);
subject_radii_mm = NaN(n_roi, n_subjects);

for subject_idx = 1:n_subjects
    for roi_idx = 1:n_roi
        % Draw independent x, y, and z translations from a Gaussian
        % distribution, restricted to within two standard deviations.
        translation = zeros(1, 3);

        for dimension_idx = 1:3
            translation(dimension_idx) = local_draw_truncated_zero_mean( ...
                sigma_mm(roi_idx), ...
                random_stream);
        end

        translated_center = centers(roi_idx, :) + translation;

        % Snap the translated center to the nearest brain node.
        squared_distance_to_center = sum( ...
            (p - translated_center).^2, ...
            2);

        [~, nearest_node_idx] = min(squared_distance_to_center);
        mesh_center = p(nearest_node_idx, :);

        % Draw a positive subject-specific radius, limited to no more than
        % two standard deviations above the requested mean radius.
        subject_radius = local_draw_radius( ...
            radius_mm(roi_idx), ...
            radius_std_mm(roi_idx), ...
            random_stream);

        % Select all nodes within the subject-specific spherical ROI.
        squared_distance_from_mesh_center = sum( ...
            (p - mesh_center).^2, ...
            2);

        roi_nodes = ...
            squared_distance_from_mesh_center <= subject_radius^2;

        subj_prob(subject_idx, roi_idx, :) = ...
            reshape(double(roi_nodes), 1, 1, n_nodes);

        center_translations(roi_idx, :, subject_idx) = translation;
        translated_centers(roi_idx, :, subject_idx) = translated_center;
        mesh_centers(roi_idx, :, subject_idx) = mesh_center;
        mesh_center_indices(roi_idx, subject_idx) = nearest_node_idx;
        subject_radii_mm(roi_idx, subject_idx) = subject_radius;
    end
end

% Preserve the nROI-by-nNodes shape when only one ROI is supplied.
avg_prob = reshape( ...
    mean(subj_prob, 1), ...
    n_roi, n_nodes);

% -------------------------------------------------------------------------
% Package diagnostics
% -------------------------------------------------------------------------
roi_info = struct();

roi_info.centers_mm = centers;
roi_info.center_translations_mm = center_translations;
roi_info.translated_centers_mm = translated_centers;
roi_info.mesh_centers_mm = mesh_centers;
roi_info.mesh_center_indices = mesh_center_indices;

roi_info.radius_mm = radius_mm;
roi_info.radius_std_mm = radius_std_mm;
roi_info.subject_radii_mm = subject_radii_mm;
roi_info.sigma_mm = sigma_mm;

roi_info.opts = opts;
end


function centers = local_resolve_centers(roi_def)
% Resolve supported ROI center field names.

if isfield(roi_def, 'centers_mm') && ~isempty(roi_def.centers_mm)
    centers = roi_def.centers_mm;

elseif isfield(roi_def, 'center_mm') && ~isempty(roi_def.center_mm)
    centers = roi_def.center_mm;

elseif isfield(roi_def, 'centers') && ~isempty(roi_def.centers)
    centers = roi_def.centers;

elseif isfield(roi_def, 'center') && ~isempty(roi_def.center)
    centers = roi_def.center;

else
    error([ ...
        'roi_def must contain center, centers, ', ...
        'center_mm, or centers_mm.']);
end

if isvector(centers) && numel(centers) == 3
    centers = reshape(centers, 1, 3);
end

if ~isnumeric(centers) || size(centers, 2) ~= 3 || ...
        any(~isfinite(centers(:)))
    error('ROI centers must be finite numeric coordinates [nROI x 3].');
end
end


function values = local_get_roi_parameter( ...
    roi_def, field_name, n_roi, default_value)
% Resolve a scalar or per-ROI parameter from roi_def.

if isfield(roi_def, field_name) && ~isempty(roi_def.(field_name))
    values = roi_def.(field_name);
else
    values = default_value;
end

if isscalar(values)
    values = repmat(values, n_roi, 1);
else
    values = values(:);

    if numel(values) ~= n_roi
        error( ...
            'roi_def.%s must be scalar or contain one value per ROI (%d).', ...
            field_name, ...
            n_roi);
    end
end
end


function value = local_draw_truncated_zero_mean(sigma, random_stream)
% Draw N(0, sigma) while restricting the result to +/-2 sigma.

if sigma == 0
    value = 0;
    return
end

limit = 2 * sigma;
value = Inf;

while abs(value) > limit
    value = sigma * local_randn(random_stream);
end
end


function radius = local_draw_radius( ...
    mean_radius, radius_std, random_stream)
% Draw a positive radius no greater than mean + 2 standard deviations.

if radius_std == 0
    radius = mean_radius;
    return
end

upper_limit = mean_radius + 2 * radius_std;
radius = -Inf;

while radius <= 0 || radius >= upper_limit
    radius = mean_radius + ...
        radius_std * local_randn(random_stream);
end
end


function value = local_randn(random_stream)
% Draw one standard-normal value from the selected random stream.

if isempty(random_stream)
    value = randn();
else
    value = randn(random_stream);
end
end