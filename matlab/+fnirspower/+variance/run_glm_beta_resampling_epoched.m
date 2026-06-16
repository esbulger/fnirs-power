function OUT = run_glm_beta_resampling_epoched(E, iter_blocks, n_iter, bad_channels, concat_baseline_sec)
%RUN_GLM_BETA_RESAMPLING_EPOCHED  Resample GLM betas from pre-epoched Hb.
%
%  OUT = RUN_GLM_BETA_RESAMPLING_EPOCHED(E, iter_blocks, n_iter, bad_channels, concat_baseline_sec)
%
%  Purpose
%  -------
%  Estimates within-subject variability of GLM beta coefficients by repeatedly
%  sampling a subset of block epochs (with replacement), concatenating them
%  into a single time series, building a block-wise design matrix, and fitting
%  an AR-IRLS GLM per channel. This mirrors the main GLM path (trend,
%  intercept, canonical HRF, optional aux regressor) to ensure consistency.
%
%  Inputs
%  ------
%  E : struct of pre-epoched data with fields
%      .hbo (B×T×nCh)  : HbO epochs (blocks × samples × channels)
%      .hbr (B×T×nCh)  : HbR epochs
%      .aux (B×T)      : auxiliary epochs (proxy or real)
%      .aux_real (opt) : if present, real aux epochs take precedence
%      .Fs             : sampling rate (Hz)
%      .T              : samples per epoch
%      .baseline_sec   : seconds of pre-stim baseline per epoch
%      .block_sec      : task duration per epoch (seconds)
%  iter_blocks : 1×nB vector of block counts to sample each iteration (e.g., [4 8 12 16])
%  n_iter      : number of resampling iterations (default 100)
%  bad_channels: optional logical/index vector for channels to set NaN
%  concat_baseline_sec : baseline length (sec) used for per-block alignment
%                        during concatenation (default = E.baseline_sec)
%
%  Outputs
%  -------
%  OUT.beta_hbo : n_iter × nB × nCh matrix of HbO betas (3rd regressor: HRF)
%  OUT.beta_hbr : n_iter × nB × nCh matrix of HbR betas
%  OUT.Fs       : sampling rate (copied from E)
%  OUT.iter_blocks, OUT.nCh : metadata
%
%  Key implementation details
%  --------------------------
%  - Concatenation uses REQUIRED per-block baseline alignment over the first
%    concat_baseline_sec seconds to avoid DC jumps at boundaries.
%  - A short half-cosine edge taper (~0.25 s) further reduces edge ripple
%    without altering within-block shapes.
%  - The design matrix is built via variance.build_concat_design, which
%    internally calls glmproc.make_design to match the main pipeline
%    (trend, intercept, canonical HRF, aux regressor).
%  - GLM is solved with glm_ar_irls (prewhitening + IRLS), Pmax ≈ 4 s.

arguments
  E struct
  iter_blocks (1,:) double
  n_iter (1,1) double = 100
  bad_channels = []
  concat_baseline_sec (1,1) double = E.baseline_sec
end

% Basic dimensions
Fs = E.Fs;                 % sampling rate (Hz)
T  = E.T;                  % samples per epoch
nCh = size(E.hbo,3);       % number of channels
if isempty(bad_channels)
  bad = [];                % no masking by default
else
  bad = bad_channels;
end

% Allocate outputs: (iterations × block-counts × channels)
beta_hbo = NaN(n_iter, numel(iter_blocks), nCh);
beta_hbr = NaN(n_iter, numel(iter_blocks), nCh);

% GLM settings: prewhitening order parameter (≈ 4 seconds)
Pmax  = ceil(4*Fs);
BL    = round(concat_baseline_sec * Fs);   % baseline samples used for alignment
taper = round(0.25 * Fs);                  % ~0.25 s edge taper (set 0 to disable)

for it = 1:n_iter
  for b = 1:numel(iter_blocks)
    B = iter_blocks(b);                    % number of blocks to sample

    % Sample block indices WITH replacement to capture variability
    idx = randi(size(E.hbo,1), B, 1);

    % === Concatenate HbO/HbR and aux epochs ===
    % Baseline alignment is mandatory to reduce boundary jumps; taper smooths edges.
    Yhbo = fnirspower.variance.concat_blocks(E.hbo(idx,:,:), BL, struct('taper_samples',taper));
    Yhbr = fnirspower.variance.concat_blocks(E.hbr(idx,:,:), BL, struct('taper_samples',taper));

    % Prefer real auxiliary epochs if provided; fall back to proxy otherwise
    if isfield(E,'aux_real')
      aux_blocks = reshape(E.aux_real(idx,:), [B T 1]);
    else
      aux_blocks = reshape(E.aux(idx,:), [B T 1]);
    end
    aux  = fnirspower.variance.concat_blocks(aux_blocks, BL, struct('taper_samples',taper));

    % === Build design matrix for the concatenated sequence ===
    % Reuses glmproc.make_design via the variance.build_concat_design wrapper
    X = fnirspower.variance.build_concat_design(Fs, T, B, E.baseline_sec, E.block_sec, aux);

    % === GLM fits (per channel) with AR-IRLS prewhitening ===
    % glm_ar_irls returns [beta; stderr; ...]; 3rd row corresponds to HRF regressor beta
    Bhbo = glm_ar_irls(Yhbo, X, Pmax); bhbo = Bhbo(3,:)';
    Bhbr = glm_ar_irls(Yhbr, X, Pmax); bhbr = Bhbr(3,:)';

    % Mask bad channels if requested (set to NaN)
    if ~isempty(bad)
      bhbo(bad) = NaN; bhbr(bad) = NaN;
    end

    % Store into iteration × block-count × channel arrays
    beta_hbo(it, b, :) = bhbo;
    beta_hbr(it, b, :) = bhbr;
  end
end

OUT = struct('beta_hbo',beta_hbo, 'beta_hbr',beta_hbr, ...
             'Fs',Fs, 'iter_blocks',iter_blocks, 'nCh',nCh);
end
