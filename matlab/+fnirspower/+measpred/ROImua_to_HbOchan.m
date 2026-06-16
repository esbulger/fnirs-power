function [HbO_subj_full, HbR_subj_full, muA760_full, muA850_full, amplitude_scale] = ROImua_to_HbOchan( ...
        p, e, brain_nodes_idx, d_mu_a_est, hb_ratio, varargin)
%FNIRSPOWER.MEASPRED.ROIMUA_TO_HBOCHAN
% Convert a subject-specific ROI map to relative Hb and absorption fields.
%
% roi_amp_std specifies the fractional standard deviation of subject amplitude
% around unit amplitude. For example, roi_amp_std = 0.3 generates:
%
%   amplitude_scale = 1 + 0.3 * randn()
%
% The same realized amplitude is applied to the relative Hb fields and to
% the wavelength-specific absorption fields used for forward projection.
%
% Outputs
% -------
% HbO_subj_full :
%     Relative HbO field over all mesh nodes [nNodes x 1].
%
% HbR_subj_full :
%     Relative HbR field over all mesh nodes [nNodes x 1].
%
% muA760_full :
%     Absorption-change field at 760 nm over all mesh nodes [nNodes x 1].
%
% muA850_full :
%     Absorption-change field at 850 nm over all mesh nodes [nNodes x 1].
%
% amplitude_scale :
%     Realized subject-level amplitude multiplier.

ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip, 'p', @isnumeric);
addRequired(ip, 'e', @isnumeric);
addRequired(ip, 'brain_nodes_idx', @isnumeric);
addRequired(ip, 'd_mu_a_est', ...
    @(x) isnumeric(x) && numel(x) == 2);
addRequired(ip, 'hb_ratio', ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(ip, 'roi_amp_std', 0.3, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(ip, 'roi_indices', [], @isnumeric);
addParameter(ip, 'roi_prob', [], @isnumeric);
addParameter(ip, 'plot_mesh', false, ...
    @(x) isnumeric(x) || islogical(x));
addParameter(ip, 'regionid', 2, ...
    @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'plot_alpha', 0.35, ...
    @(x) isnumeric(x) && isscalar(x));

parse( ...
    ip, ...
    p, e, brain_nodes_idx, d_mu_a_est, hb_ratio, ...
    varargin{:});

opts = ip.Results;

% Use the parsed and validated required inputs.
p = opts.p;
e = opts.e;
brain_nodes_idx = opts.brain_nodes_idx(:);
d_mu_a_est = opts.d_mu_a_est(:);
hb_ratio = opts.hb_ratio;

% Accept roi_prob as an alias for roi_indices.
roi_indices = opts.roi_indices;

if isempty(roi_indices) && ~isempty(opts.roi_prob)
    roi_indices = opts.roi_prob;
end

if isempty(roi_indices)
    error( ...
        'Provide roi_indices or roi_prob as a vector over brain nodes.');
end

roi_indices = roi_indices(:);

if numel(roi_indices) ~= numel(brain_nodes_idx)
    error( ...
        ['roi_indices must contain one value per brain node. ', ...
         'Expected %d values but received %d.'], ...
        numel(brain_nodes_idx), ...
        numel(roi_indices));
end

% Normalize the spatial template to unit peak so that amplitude_scale
% controls the peak response.
template = roi_indices;
template = template / max(eps, max(abs(template)));

% Generate one amplitude realization for this subject and ROI.
amplitude_scale = 1 + opts.roi_amp_std * randn();

% Construct relative voxelwise Hb fields.
HbO_subj = amplitude_scale * template;
HbR_subj = -(amplitude_scale / hb_ratio) * template;

HbO_subj_full = zeros(size(p, 1), 1);
HbR_subj_full = zeros(size(p, 1), 1);

HbO_subj_full(brain_nodes_idx) = HbO_subj;
HbR_subj_full(brain_nodes_idx) = HbR_subj;

% Construct wavelength-specific absorption fields using the same spatial
% template and subject-level amplitude realization.
muA760_full = zeros(size(p, 1), 1);
muA850_full = zeros(size(p, 1), 1);

muA760_full(brain_nodes_idx) = ...
    amplitude_scale * template * d_mu_a_est(1);

muA850_full(brain_nodes_idx) = ...
    amplitude_scale * template * d_mu_a_est(2);

if opts.plot_mesh
    tet = e(e(:, 5) == opts.regionid, :);

    figure;
    hold on;

    iso2mesh_plotmesh( ...
        [p HbO_subj_full], ...
        tet, ...
        'EdgeAlpha', 0.1, ...
        'FaceAlpha', opts.plot_alpha);

    title('Relative HbO ROI field');
    axis equal;
    axis off;
    colorbar;
end
end