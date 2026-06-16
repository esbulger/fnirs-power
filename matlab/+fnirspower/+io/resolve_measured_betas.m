function [meas_HbO, meas_HbR, have_group, all_HbO, all_HbR] = ...
    resolve_measured_betas(B)
%FNIRSPOWER.IO.RESOLVE_MEASURED_BETAS
% Resolve measured channel-level HbO and HbR beta estimates.
%
% Supported inputs
% ----------------
% Unmasked GLM beta file:
%
%   B.GLM_Hb.all_subj_beta_hbo
%   B.GLM_Hb.all_subj_beta_hbr
%   B.GLM_Hb.mean_hbo
%   B.GLM_Hb.mean_hbr
%
% Masked group-mean beta file:
%
%   B.GLM_Hb_Masked.mean_hbo
%   B.GLM_Hb_Masked.mean_hbr
%
% The GLM_Hb or GLM_Hb_Masked structure may also be passed directly.
%
% Outputs
% -------
% meas_HbO, meas_HbR :
%     Group-mean channel beta vectors [1 x nChannels].
%
% have_group :
%     True when per-subject beta arrays are available.
%
% all_HbO, all_HbR :
%     Per-subject channel beta arrays [nSubjects x nChannels].
%     These are empty when a masked group-mean file is supplied.

if ~isstruct(B) || ~isscalar(B)
    error('B must be a scalar structure.');
end

have_group = false;
all_HbO = [];
all_HbR = [];

% Resolve the structure stored in the loaded MAT file.
if isfield(B, 'GLM_Hb') && isfield(B, 'GLM_Hb_Masked')
    error('Input contains both GLM_Hb and GLM_Hb_Masked structures.');

elseif isfield(B, 'GLM_Hb')
    if ~isstruct(B.GLM_Hb) || ~isscalar(B.GLM_Hb)
        error('B.GLM_Hb must be a scalar structure.');
    end

    G = B.GLM_Hb;

elseif isfield(B, 'GLM_Hb_Masked')
    if ~isstruct(B.GLM_Hb_Masked) || ...
            ~isscalar(B.GLM_Hb_Masked)
        error('B.GLM_Hb_Masked must be a scalar structure.');
    end

    G = B.GLM_Hb_Masked;

else
    % Permit passing GLM_Hb or GLM_Hb_Masked directly.
    G = B;
end

has_hbo_subjects = isfield(G, 'all_subj_beta_hbo');
has_hbr_subjects = isfield(G, 'all_subj_beta_hbr');

if xor(has_hbo_subjects, has_hbr_subjects)
    error([ ...
        'Per-subject HbO and HbR fields must either both be present ', ...
        'or both be absent.']);
end

has_subject_betas = has_hbo_subjects && has_hbr_subjects;

has_hbo_mean = isfield(G, 'mean_hbo') && ~isempty(G.mean_hbo);
has_hbr_mean = isfield(G, 'mean_hbr') && ~isempty(G.mean_hbr);

if xor(has_hbo_mean, has_hbr_mean)
    error([ ...
        'mean_hbo and mean_hbr must either both be present ', ...
        'or both be absent.']);
end

has_saved_means = has_hbo_mean && has_hbr_mean;

% Resolve and validate per-subject beta arrays when available.
if has_subject_betas
    all_HbO = G.all_subj_beta_hbo;
    all_HbR = G.all_subj_beta_hbr;

    if ~isnumeric(all_HbO) || ~isnumeric(all_HbR) || ...
            ~ismatrix(all_HbO) || ~ismatrix(all_HbR)
        error( ...
            'Per-subject HbO and HbR beta arrays must be numeric matrices.');
    end

    if isempty(all_HbO) || isempty(all_HbR)
        error('Per-subject HbO and HbR beta arrays cannot be empty.');
    end

    if ~isequal(size(all_HbO), size(all_HbR))
        error([ ...
            'all_subj_beta_hbo and all_subj_beta_hbr must have ', ...
            'the same dimensions.']);
    end

    have_group = true;
end

% Use saved means when present. Otherwise calculate them from subjects.
if has_saved_means
    meas_HbO = G.mean_hbo(:)';
    meas_HbR = G.mean_hbr(:)';

elseif has_subject_betas
    meas_HbO = mean(all_HbO, 1, 'omitnan');
    meas_HbR = mean(all_HbR, 1, 'omitnan');

else
    error([ ...
        'Unrecognized GLM beta structure. Expected mean_hbo and ', ...
        'mean_hbr, optionally with all_subj_beta_hbo and ', ...
        'all_subj_beta_hbr.']);
end

if ~isnumeric(meas_HbO) || ~isnumeric(meas_HbR)
    error('mean_hbo and mean_hbr must be numeric.');
end

if numel(meas_HbO) ~= numel(meas_HbR)
    error('The resolved HbO and HbR mean vectors must have equal lengths.');
end

if have_group && numel(meas_HbO) ~= size(all_HbO, 2)
    error([ ...
        'The number of group-mean channels does not match the ', ...
        'per-subject beta arrays.']);
end
end