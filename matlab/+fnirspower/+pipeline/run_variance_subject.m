function OUT = run_variance_subject(snirf_path, nirs_model, baseline_sec, block_sec, iter_blocks, n_iter, concat_baseline_sec)
%RUN_VARIANCE_SUBJECT Run subject-level block-resampling analysis for GLM beta variability.
%
%  OUT = fnirspower.pipeline.run_variance_subject( ...
%      snirf_path, nirs_model, baseline_sec, block_sec, iter_blocks, n_iter, concat_baseline_sec)
%
%  Description
%  -----------
%  This function performs the subject-level portion of the variance analysis
%  workflow. It:
%
%    1) Loads and preprocesses a single SNIRF recording.
%    2) Epochs the recording into blocks using the specified baseline and
%       block durations.
%    3) Repeatedly resamples blocks for each requested block-count condition.
%    4) Concatenates the resampled blocks into synthetic time series.
%    5) Refits GLMs to estimate the variability of channel-level HbO/HbR
%       beta values across resampling iterations.
%
%  This function is typically called internally by
%  fnirspower.pipeline.run_variance, but it can also be used directly for
%  subject-level variance analysis.
%
%  Required inputs
%  ---------------
%  snirf_path :
%      Path to a single subject/session SNIRF file.
%
%  nirs_model:
%      Forward-model/preprocessing struct passed to
%      fnirspower.nirsproc.preprocess. This should contain the information
%      needed for the standard preprocessing workflow.
%
%  baseline_sec :
%      Baseline duration in seconds used during block extraction.
%
%  block_sec :
%      Duration of each task block in seconds.
%
%  iter_blocks :
%      Vector of block counts to evaluate during resampling, e.g.
%      [4 8 12 16 20 24].
%
%  n_iter :
%      Number of resampling iterations to run for each block-count condition.
%
%  concat_baseline_sec :
%      Baseline duration in seconds used during concatenation of resampled
%      blocks.
%
%  Returns
%  -------
%  OUT : struct with fields
%
%    OUT.beta_hbo :
%        HbO beta estimates across iterations and block-count conditions:
%          [nIter x nBlockConds x nCh]
%
%    OUT.beta_hbr :
%        HbR beta estimates across iterations and block-count conditions:
%          [nIter x nBlockConds x nCh]
%
%    OUT.Fs :
%        Sampling frequency of the preprocessed data.
%
%    OUT.nCh :
%        Number of channels included in the analysis.
%
%    OUT.iter_blocks :
%        Vector of block counts evaluated.
%
%    OUT.bad_channels :
%        Logical or index vector of bad channels identified during
%        preprocessing.
%
%  Example
%  -------
%  OUT = fnirspower.pipeline.run_variance_subject( ...
%      snirf_path, model, 2.5, 15, [4 8 12 16 20 24], 100, 2.5);
%
%  Notes
%  -----
%  - Preprocessing is performed once at the start of the function using
%    fnirspower.nirsproc.preprocess.
%  - Block resampling and GLM refitting are delegated to
%    fnirspower.variance.run_glm_beta_resampling.
%  - This function does not compute group-level summaries; it returns only
%    the subject-level resampling results needed by
%    fnirspower.pipeline.run_variance.
%
%  See also
%  --------
%  fnirspower.pipeline.run_variance
%  fnirspower.nirsproc.preprocess
%  fnirspower.variance.run_glm_beta_resampling

arguments
  snirf_path (1,:) char
  nirs_model struct
  baseline_sec (1,1) double
  block_sec (1,1) double
  iter_blocks (1,:) double
  n_iter (1,1) double
  concat_baseline_sec (1,1) double
end

% 1) preprocess (same as GLM pipeline)
S = fnirspower.nirsproc.preprocess(snirf_path, nirs_model);

% 2) resampling GLM on epoched/concatenated blocks
fprintf("Initializing resampling and GLM for subject from %s", snirf_path)
R = fnirspower.variance.run_glm_beta_resampling( ...
      S, baseline_sec, block_sec, iter_blocks, n_iter, S.bad_channels, concat_baseline_sec);

OUT = struct();
OUT.beta_hbo = R.beta_hbo;     % n_iter x nB x nCh
OUT.beta_hbr = R.beta_hbr;     % n_iter x nB x nCh
OUT.Fs = R.Fs;
OUT.nCh = R.nCh;
OUT.iter_blocks = R.iter_blocks;
OUT.bad_channels = S.bad_channels;
end