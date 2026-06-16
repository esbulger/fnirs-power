function G = run_group_glm(subject_ids, snirf_paths, mesh_path, nirsmodel_path, layout_path, opts)
%RUN_GROUP_GLM Run subject-level GLM and source reconstruction across multiple SNIRF recordings.
%
%  G = fnirspower.pipeline.run_group_glm( ...
%      subject_ids, paths, mesh_path, nirsmodel_path, layout_path, ...)
%
%  Description
%  -----------
%  This function performs group-level orchestration of the subject GLM
%  workflow. It:
%
%    1) Validates one SNIRF file per subject from paths.snirf_files.
%    2) Loads the shared mesh and forward-model data used for reconstruction.
%    3) Builds the HbO/HbR Jacobians from the wavelength-specific Jacobians.
%    4) Runs fnirspower.pipeline.run_subject_glm once per subject.
%    5) Collects channel-level betas, source-space reconstructions, and
%       reconstructed spectral absorption maps across subjects.
%    6) Computes group-average spectral absorption maps and one-sample
%       t-tests against zero.
%
%  Required inputs
%  ---------------
%  subject_ids :
%      Numeric vector of subject IDs. The order of subject_ids must match
%      the order of paths.snirf_files.
%
%  snirf_paths:
%      Struct containing:
%
%        snirf_paths.snirf_files :
%            Cell array or string array of full SNIRF file paths, one per
%            subject, in the same order as subject_ids.
%
%  mesh_path :
%      Path to the tetrahedral head mesh MAT file. Used to recover brain-node
%      indices from either:
%        - genmesh.ele
%      or:
%        - e
%
%  nirsmodel_path :
%      Path to the forward-model MAT file. Expected to contain:
%
%        J_760.complete :
%            Full 760 nm Jacobian
%
%        J_850.complete :
%            Full 850 nm Jacobian
%
%        nirs_mesh :
%            NIRFAST/NIRFASTer mesh structure used during preprocessing and
%            downstream reconstruction
%
%  layout_path:
%      Optional layout MAT file. If supplied and non-empty, the raw loaded
%      layout struct is attached to the output as G.layout.
%
%  Name-value options
%  ------------------
%  'hrf_seconds' :
%      Duration of the modeled hemodynamic response function in seconds.
%      Passed through to fnirspower.pipeline.run_subject_glm.
%      Default: 32
%
%  'alpha' :
%      Regularization parameter used during source reconstruction.
%      Passed through to fnirspower.pipeline.run_subject_glm.
%      Default: 0.01
%
%  'n_pca' :
%      Number of PCA nuisance regressors retained during subject-level
%      preprocessing / design construction.
%      Default: 1
%
%  'r2max' :
%      Maximum R^2 threshold used during PCA nuisance regressor selection.
%      Default: 2.5
%
%  'Pmax' :
%      Optional maximum AR model order passed through to the subject-level
%      GLM fitting routine.
%      Default: []
%
%  Returns
%  -------
%  G : struct with fields
%
%    G.subjects :
%        Subject IDs included in the group run.
%
%    G.snirf_files :
%        SNIRF file paths used for each subject.
%
%    G.beta :
%        Struct containing channel-level beta arrays with one row per subject:
%
%          .hbo   :
%              HbO betas [nSubj x nCh]
%
%          .hbr   :
%              HbR betas [nSubj x nCh]
%
%          .od760 :
%              760 nm ΔOD betas [nSubj x nCh]
%
%          .od850 :
%              850 nm ΔOD betas [nSubj x nCh]
%
%    G.recon :
%        Struct containing source-space hemoglobin reconstructions:
%
%          .hbo :
%              HbO reconstructions [nSubj x nV]
%
%          .hbr :
%              HbR reconstructions [nSubj x nV]
%
%    G.mua :
%        Struct containing reconstructed spectral absorption changes on
%        brain nodes:
%
%          .lambda760 :
%              Δμa at 760 nm [nSubj x nBrain]
%
%          .lambda850 :
%              Δμa at 850 nm [nSubj x nBrain]
%
%    G.group :
%        Struct containing group-level spectral absorption summaries:
%
%          .mean_mua760 :
%              Mean Δμa760 across subjects
%
%          .mean_mua850 :
%              Mean Δμa850 across subjects
%
%          .pvals_mua760 :
%              One-sample t-test p-values for Δμa760 against zero
%
%          .pvals_mua850 :
%              One-sample t-test p-values for Δμa850 against zero
%
%  Example
%  -------
%  snirf_paths = struct();
%  snirf_paths.snirf_files = {
%      '/abs/path/sub101_run.snirf'
%      '/abs/path/sub102_run.snirf'
%      '/abs/path/sub103_run.snirf'
%  };
%
%  G = fnirspower.pipeline.run_group_glm( ...
%      [101 102 103], snirf_paths, mesh_path, nirsmodel_path, layout_path, ...
%      'hrf_seconds', 32, ...
%      'alpha', 0.01, ...
%      'n_pca', 1, ...
%      'r2max', 2.5, ...
%      'Pmax', []);
%
%  Notes
%  -----
%  - subject_ids and paths.snirf_files must be aligned in the same order.
%  - This function assumes a shared forward model across all subjects in the
%    group run.
%  - Group-level statistics are currently computed only for reconstructed
%    spectral absorption maps (Δμa760 and Δμa850), not for channel-level
%    HbO/HbR betas.
%  - Subject-level preprocessing, GLM fitting, and reconstruction are
%    delegated to fnirspower.pipeline.run_subject_glm.
%
%  See also
%  --------
%  fnirspower.pipeline.run_subject_glm
%  fnirspower.recon.invert_woodbury

arguments
    subject_ids (1,:) double
    snirf_paths struct
    mesh_path (1,:) char
    nirsmodel_path (1,:) char
    layout_path (1,:) char = ''
    opts.hrf_seconds double = 32
    opts.alpha double = 0.01
    opts.n_pca double = 1
    opts.r2max double = 2.5
    opts.Pmax double = []
end

if ~isfield(snirf_paths, 'snirf_files') || isempty(snirf_paths.snirf_files)
    error('paths.snirf_files is required and must match subject_ids.');
end

snirf_files = snirf_paths.snirf_files;
if isstring(snirf_files)
    snirf_files = cellstr(snirf_files(:));
end
if ~iscell(snirf_files)
    error('paths.snirf_files must be a cell array or string array of full SNIRF paths.');
end
if numel(snirf_files) ~= numel(subject_ids)
    error('paths.snirf_files must have the same length/order as subject_ids.');
end

M = load(mesh_path);
if isfield(M,'genmesh')
    e = M.genmesh.ele;
else
    e = M.e;
end


nirs_model = load(nirsmodel_path);

J760_full = abs(nirs_model.J_760.complete);
J850_full = abs(nirs_model.J_850.complete);
J_HbO = J760_full * 58.6/(58.6+105.8) + J850_full * 105.8/(58.6+105.8);
J_HbR = J760_full * 154.8/(69.1+154.8) + J850_full *  69.1/(69.1+154.8);

brain_tet = e(e(:,5)==2,:);
brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);
nirs_model.epsilon_mm_uM = 2.303e-6 * [58.6 154.8; 105.8 69.1];
nirs_model.brain_nodes_idx = brain_nodes_idx;
J = struct('J760_full',J760_full,'J850_full',J850_full,'J_HbO',J_HbO,'J_HbR',J_HbR);

G = struct();
G.subjects = subject_ids;
G.snirf_files = snirf_files;

for i = 1:numel(subject_ids)
    sid = subject_ids(i);
    snirf_path = char(snirf_files{i});
    assert(exist(snirf_path,'file')==2, 'SNIRF file not found for subject %d: %s', sid, snirf_path);

    fprintf('Processing subject %d...\n', sid);
    S_out = fnirspower.pipeline.run_subject_glm( ...
        snirf_path, nirs_model, J, ...
        'hrf_seconds', opts.hrf_seconds, ...
        'alpha', opts.alpha, ...
        'n_pca', opts.n_pca, ...
        'r2max', opts.r2max, ...
        'Pmax', opts.Pmax);

    if i==1
        nCh = numel(S_out.beta.hbo);
        nV  = numel(S_out.recon.hbo);
        nBrain = numel(S_out.mua.lambda850);

        G.beta.hbo   = NaN(numel(subject_ids), nCh);
        G.beta.hbr   = NaN(numel(subject_ids), nCh);
        G.beta.od760 = NaN(numel(subject_ids), nCh);
        G.beta.od850 = NaN(numel(subject_ids), nCh);

        G.recon.hbo  = NaN(numel(subject_ids), nV);
        G.recon.hbr  = NaN(numel(subject_ids), nV);

        G.mua.lambda760 = NaN(numel(subject_ids), nBrain);
        G.mua.lambda850 = NaN(numel(subject_ids), nBrain);
    end

    G.beta.hbo(i,:)      = S_out.beta.hbo.';
    G.beta.hbr(i,:)      = S_out.beta.hbr.';
    G.beta.od760(i,:)    = S_out.beta.od760.';
    G.beta.od850(i,:)    = S_out.beta.od850.';
    G.recon.hbo(i,:)     = S_out.recon.hbo.';
    G.recon.hbr(i,:)     = S_out.recon.hbr.';
    G.mua.lambda760(i,:) = S_out.mua.lambda760.';
    G.mua.lambda850(i,:) = S_out.mua.lambda850.';
end

G.group.mean_mua760 = mean(G.mua.lambda760, 1, 'omitnan');
G.group.mean_mua850 = mean(G.mua.lambda850, 1, 'omitnan');
[~,G.group.pvals_mua760] = ttest(G.mua.lambda760, 0);
[~,G.group.pvals_mua850] = ttest(G.mua.lambda850, 0);

end