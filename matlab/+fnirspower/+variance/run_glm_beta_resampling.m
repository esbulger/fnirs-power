function OUT = run_glm_beta_resampling(S, baseline_sec, block_sec, iter_blocks, n_iter, bad_channels, concat_baseline_sec)
%RUN_GLM_BETA_RESAMPLING  Resample GLM betas from raw SNIRF via preprocessing.
%
%  OUT = RUN_GLM_BETA_RESAMPLING(S, baseline_sec, block_sec, iter_blocks, n_iter, bad_channels, concat_baseline_sec)
%
%  Purpose
%  -------
%  Estimates within-subject variability of GLM beta coefficients by repeatedly
%  sampling a subset of block epochs (with replacement) from a single subject,
%  concatenating them into a continuous time series, building a block-consistent
%  design matrix, and fitting an AR-IRLS GLM per channel. This mirrors the main
%  GLM path (trend, intercept, canonical HRF, optional aux) to ensure consistency.
%
%  Inputs
%  ------
%  S : struct from nirsproc.preprocess for one subject, containing fields:
%      .HBO, .HBD, .dOD_hp, .Fs, .time, .stimvec, .nCh, and optionally aux
%  baseline_sec : seconds of pre-stim baseline per epoch (used for epoching)
%  block_sec    : seconds of task duration per epoch
%  iter_blocks  : 1×nB vector of block counts to sample each iteration (e.g., [4 8 12 16])
%  n_iter       : number of resampling iterations (default 100)
%  bad_channels : optional logical/index vector for channels to set NaN (default S.bad_channels)
%  concat_baseline_sec : baseline length (sec) used for per-block alignment
%                        during concatenation (default = baseline_sec)
%
%  Outputs
%  -------
%  OUT.beta_hbo : n_iter × nB × nCh matrix of HbO betas (HRF regressor)
%  OUT.beta_hbr : n_iter × nB × nCh matrix of HbR betas
%  OUT.Fs       : sampling rate
%  OUT.iter_blocks, OUT.nCh : metadata

arguments
    S struct
    baseline_sec (1,1) double
    block_sec (1,1) double
    iter_blocks (1,:) double
    n_iter (1,1) double = 100
    bad_channels = []
    concat_baseline_sec (1,1) double = baseline_sec
end

% Epoch HbO and HbR. Use the measured accelerometer regressor when one was
% recovered during preprocessing; otherwise, hb_epochs constructs its
% optical-density-based proxy.
if isfield(S, 'has_real_acc') && S.has_real_acc
    if ~isfield(S, 'acc_reg') || isempty(S.acc_reg)
        error( ...
            'RUN_GLM_BETA_RESAMPLING:MissingAccelerometer', ...
            ['S.has_real_acc is true, but S.acc_reg is missing or empty. ', ...
             'Check the preprocessing output.']);
    end

    E = fnirspower.nirsproc.hb_epochs( ...
        S, ...
        baseline_sec, ...
        block_sec, ...
        S.acc_reg);
else
    E = fnirspower.nirsproc.hb_epochs( ...
        S, ...
        baseline_sec, ...
        block_sec);
end

% Basic dimensions
Fs = E.Fs;                   % sampling rate (Hz)
T  = E.T;                    % samples per epoch
nCh = size(E.hbo,3);         % number of channels
if isempty(bad_channels)
    bad = S.bad_channels;      % default to subject-specific mask
else
    bad = bad_channels;
end

% Allocate outputs: (iterations × block-counts × channels)
beta_hbo = NaN(n_iter, numel(iter_blocks), nCh);
beta_hbr = NaN(n_iter, numel(iter_blocks), nCh);

% GLM settings
Pmax  = ceil(4*Fs);                  % AR prewhitening order ~ 4 seconds
BL    = round(concat_baseline_sec * Fs);  % REQUIRED baseline alignment window (samples)
taper = round(0.25 * Fs);           % ~0.25 s edge taper (set 0 to disable)

for it = 1:n_iter
    fprintf("Iteration: %i / %i", it, n_iter)
    for b = 1:numel(iter_blocks)
        B = iter_blocks(b);             % number of blocks to sample
        
        % Sample block indices WITH replacement to capture variability
        idx = randi(size(E.hbo,1), B, 1);
        
        % === Concatenate HbO/HbR and aux epochs ===
        % Per-block baseline alignment (first BL samples per channel) is
        % mandatory to reduce boundary DC jumps; taper smooths edges without
        % warping within-block shapes.
        Yhbo = fnirspower.variance.concat_blocks(E.hbo(idx,:,:), BL, taper);
        Yhbr = fnirspower.variance.concat_blocks(E.hbr(idx,:,:), BL, taper);
        
        % Prefer real auxiliary epochs if provided; fall back to proxy otherwise
        if isfield(E,'aux_real')
            aux_blocks = reshape(E.aux_real(idx,:), [B T 1]);
        else
            aux_blocks = reshape(E.aux(idx,:), [B T 1]);
        end
        aux  = fnirspower.variance.concat_blocks(aux_blocks, BL, taper);
        
        % === Build design matrix for the concatenated sequence ===
        % Reuse glmproc.make_design via the variance.build_concat_design wrapper
        X = fnirspower.variance.build_concat_design(Fs, T, B, baseline_sec, block_sec, aux, 'hrf_seconds', 32, 'normalize_cols', true);
        
        % === GLM fits (per channel) with AR-IRLS prewhitening ===
        % glm_ar_irls returns rows for regressors; 3rd row corresponds to the
        % HRF regressor beta by construction of glmproc.make_design
        Bhbo = glm_ar_irls(Yhbo, X, Pmax); bhbo = Bhbo(3,:)';
        Bhbr = glm_ar_irls(Yhbr, X, Pmax); bhbr = Bhbr(3,:)';
        
        % Mask bad channels if requested (set to NaN)
        bhbo(bad) = NaN; bhbr(bad) = NaN;
        
        % Store into iteration × block-count × channel arrays
        beta_hbo(it, b, :) = bhbo;
        beta_hbr(it, b, :) = bhbr;
    end
end

OUT = struct('beta_hbo',beta_hbo, 'beta_hbr',beta_hbr, ...
    'Fs',Fs, 'iter_blocks',iter_blocks, 'nCh',nCh);
end