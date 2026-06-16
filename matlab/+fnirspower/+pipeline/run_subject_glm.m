function out = run_subject_glm(snirf_path, nirs_model, jacobians, opts)
%RUN_SUBJECT_GLM Run preprocessing, channel-level GLM fitting, and source reconstruction for one SNIRF recording.
%
%  out = fnirspower.pipeline.run_subject_glm( ...
%      snirf_path, nirs_model, jacobians, ...)
%
%  Description
%  -----------
%  This function performs the full single-subject AR-ILS GLM and reconstruction
%  workflow for one SNIRF recording. It:
%
%    1) Loads and preprocesses the SNIRF data.
%    2) Builds nuisance regressors using PCA-based components.
%    3) Constructs design matrices for HbO/HbR and ΔOD analyses.
%    4) Fits channel-level GLMs for HbO, HbR, 760 nm ΔOD, and 850 nm ΔOD.
%    5) Reconstructs source-space HbO, HbR, and spectral absorption changes.
%    6) Returns channel-level betas, reconstructed source maps, and basic
%       subject/session metadata.
%
%  Required inputs
%  ---------------
%  snirf_path :
%      Path to a subject/session SNIRF file to analyze.
%
%  nirs_model :
%      Struct describing the forward model and reconstruction metadata.
%      Expected fields include:
%
%        model.nirs_mesh :
%            NIRFAST/NIRFASTer mesh structure used during preprocessing.
%
%        model.epsilon_mm_uM :
%            Extinction coefficient matrix used for spectral conversions.
%
%        model.brain_nodes_idx :
%            Indices of brain nodes retained for source-space outputs.
%
%  jacobians :
%      Struct containing the Jacobians needed for channel-level inversion
%      and source reconstruction. Expected fields:
%
%        jacobians.J760_full :
%            Full 760 nm Jacobian.
%
%        jacobians.J850_full :
%            Full 850 nm Jacobian.
%
%        jacobians.J_HbO :
%            HbO-space Jacobian.
%
%        jacobians.J_HbR :
%            HbR-space Jacobian.
%
%  Name-value options
%  ------------------
%  'hrf_seconds' :
%      Duration of the modeled hemodynamic response function in seconds.
%      Passed into design-matrix construction. Default: 32
%
%  'alpha' :
%      Regularization parameter used during source reconstruction.
%      Passed to fnirspower.recon.invert_woodbury. Default: 0.01
%
%  'n_pca' :
%      Number of PCA nuisance regressors to retain for Hb and ΔOD design
%      matrices. Default: 1
%
%  'r2max' :
%      Maximum R^2 threshold used during PCA regressor selection.
%      Default: 2.5
%
%  'Pmax' :
%      Optional maximum AR model order passed to the GLM fitting routine.
%      Default: []
%
%  Returns
%  -------
%  out : struct with fields
%
%    out.beta :
%        Struct containing channel-level GLM beta estimates:
%
%          .hbo   :
%              HbO beta vector [nCh x 1]
%
%          .hbr   :
%              HbR beta vector [nCh x 1]
%
%          .od760 :
%              760 nm ΔOD beta vector [nCh x 1]
%
%          .od850 :
%              850 nm ΔOD beta vector [nCh x 1]
%
%    out.recon :
%        Struct containing source-space reconstructions:
%
%          .hbo :
%              Reconstructed HbO source map
%
%          .hbr :
%              Reconstructed HbR source map
%
%    out.mua :
%        Struct containing reconstructed spectral absorption changes on
%        brain nodes only:
%
%          .lambda760 :
%              Reconstructed Δμa at 760 nm [nBrain x 1]
%
%          .lambda850 :
%              Reconstructed Δμa at 850 nm [nBrain x 1]
%
%    out.info :
%        Struct containing basic subject/session analysis metadata:
%
%          .bad_channels :
%              Logical or index vector of bad channels identified during
%              preprocessing
%
%          .Fs :
%              Sampling frequency
%
%          .time :
%              Time vector used for GLM construction
%
%          .stimvec :
%              Stimulus vector used for design-matrix construction
%
%          .nCh :
%              Number of channels
%
%  Example
%  -------
%  out = fnirspower.pipeline.run_subject_glm( ...
%      snirf_path, nirs_model, jacobians, ...
%      'hrf_seconds', 32, ...
%      'alpha', 0.01, ...
%      'n_pca', 1, ...
%      'r2max', 2.5, ...
%      'Pmax', []);
%
%  Notes
%  -----
%  - This function is the main one-subject SNIRF-to-GLM pipeline entry point.
%  - Preprocessing is delegated to fnirspower.nirsproc.preprocess.
%  - PCA nuisance regressor selection is delegated to
%    fnirspower.nirsproc.pick_pca_regs.
%  - Design-matrix construction is delegated to fnirspower.glmproc.make_design.
%  - Source reconstruction uses fnirspower.recon.invert_woodbury with the
%    supplied Jacobians and reconstruction regularization parameter.
%  - HbO/HbR bad channels are explicitly set to NaN before reconstruction.
%
%  See also
%  --------
%  fnirspower.nirsproc.preprocess
%  fnirspower.nirsproc.pick_pca_regs
%  fnirspower.glmproc.make_design
%  fnirspower.glmproc.fit_glm
%  fnirspower.recon.invert_woodbury

arguments
  snirf_path (1,:) char
  nirs_model struct
  jacobians struct
  opts.hrf_seconds double = 32
  opts.alpha double = 0.01
  opts.n_pca double = 1
  opts.r2max double = 2.5
  opts.Pmax double = []
end

% 1) preprocess
S = fnirspower.nirsproc.preprocess(snirf_path, nirs_model);

% 2) PCA regs
[pca_Hb, pca_OD] = fnirspower.nirsproc.pick_pca_regs(S.HBO, S.HBD, S.dOD_hp, S.stimvec, opts.n_pca, opts.r2max);

% 3) Design matrices
X_Hb = fnirspower.glmproc.make_design(S.time, S.Fs, S.stimvec, S.acc_reg, pca_Hb, 'hrf_seconds', opts.hrf_seconds, 'normalize_cols',true);
X_OD = fnirspower.glmproc.make_design(S.time, S.Fs, S.stimvec, S.acc_reg, pca_OD, 'hrf_seconds', opts.hrf_seconds, 'normalize_cols',true);

% 4) GLM fits
fprintf('Fitting GLMs... \n')
B_hbo = fnirspower.glmproc.fit_glm(S.HBO,     X_Hb, S.Fs, 'Pmax',opts.Pmax); beta_hbo = B_hbo(3,:)';
B_hbr = fnirspower.glmproc.fit_glm(S.HBD,     X_Hb, S.Fs, 'Pmax',opts.Pmax); beta_hbr = B_hbr(3,:)';
B_760 = fnirspower.glmproc.fit_glm(S.dOD_760, X_OD, S.Fs, 'Pmax',opts.Pmax); beta_760 = B_760(3,:)';
B_850 = fnirspower.glmproc.fit_glm(S.dOD_850, X_OD, S.Fs, 'Pmax',opts.Pmax); beta_850 = B_850(3,:)';

beta_hbo(S.bad_channels) = NaN; beta_hbr(S.bad_channels) = NaN;

% 5) Recon
fprintf('Reconstructing sources... \n')
recon_hbo = fnirspower.recon.invert_woodbury(jacobians.J_HbO,      beta_hbo, S.bad_channels, opts.alpha);
recon_hbr = fnirspower.recon.invert_woodbury(jacobians.J_HbR,      beta_hbr, S.bad_channels, opts.alpha);
mu_a_760  = fnirspower.recon.invert_woodbury(jacobians.J760_full,  beta_760, S.bad_channels, opts.alpha);
mu_a_850  = fnirspower.recon.invert_woodbury(jacobians.J850_full,  beta_850, S.bad_channels, opts.alpha);

% 6) Package
out = struct();
out.beta = struct('hbo',beta_hbo,'hbr',beta_hbr,'od760',beta_760,'od850',beta_850);
out.recon = struct('hbo',recon_hbo,'hbr',recon_hbr);
out.mua  = struct('lambda760', mu_a_760(nirs_model.brain_nodes_idx), ...
                  'lambda850', mu_a_850(nirs_model.brain_nodes_idx));
out.info = struct('bad_channels',S.bad_channels,'Fs',S.Fs,'time',S.time,'stimvec',S.stimvec,'nCh',S.nCh);
end
