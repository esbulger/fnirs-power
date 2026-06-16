function out = run_absorption_estimation(muA_file, nirsmodel_path, mesh_mat, beta_file, layout_file, opts)
%RUN_ABSORPTION_ESTIMATION Estimate spectral absorption magnitude from reconstructed Δμa and measured GLM betas.
%
%  out = fnirspower.pipeline.run_absorption_estimation( ...
%      muA_file, nirsmodel_path, mesh_mat, beta_file, layout_file, ...)
%
%  Description
%  -----------
%  This function performs spectral absorption magnitude estimation by combining
%  reconstructed spectral absorption-change maps (Δμa) with measured channel-
%  level GLM betas. It:
%
%    1) Loads an existing Δμa reconstruction, or reconstructs Δμa from ΔOD
%       betas if no μa file is supplied.
%    2) Builds a spectral/chromophore-space forward model from the supplied
%       NIRFAST model and extinction coefficients.
%    3) Forward-projects each subject’s Δμa maps to the channel level.
%    4) Compares predicted channel responses to measured HbO/HbR GLM betas.
%    5) Thresholds channels based on measured values and prediction statistics.
%    6) Selects a regularization weight (lambda) by L-curve analysis.
%    7) Solves for the best-fitting spectral absorption magnitude estimate.
%    8) Produces channel-level estimated HbO/HbR maps and optional plots.
%
%  Required inputs
%  ---------------
%  muA_file :
%      Path to a MAT file containing reconstructed spectral absorption maps.
%      Expected variables:
%        - all_mu_a760   [nSubj x nBrain]
%        - all_mu_a850   [nSubj x nBrain]
%      May be passed as '' to trigger reconstruction from ΔOD betas via
%      opts.beta_od_file.
%
%  nirsmodel_path :
%      Path to a forward-model MAT file. Used to build the chromophore-space
%      spectral Jacobian through fnirspower.fwdcomp.buildPathAdjJacobian.
%
%  mesh_mat :
%      Path to the tetrahedral mesh MAT file. Used to recover brain-node
%      indices if they are not present in muA_file.
%
%  beta_file :
%      Path to a measured GLM beta MAT file. This file must contain the
%      channel-level HbO/HbR beta values required by
%      fnirspower.io.resolve_measured_betas.
%
%  layout_file :
%      Path to a layout MAT file used for optional channel-level plotting.
%      Pass '' to skip plotting.
%
%  Name-value options
%  ------------------
%  'epsilon_mm_uM' :
%      Extinction coefficient matrix used to build the spectral Jacobian.
%      Default:
%        2.303e-6 * [58.6,154.8; 105.8,69.1]
%
%  'alpha' :
%      Significance threshold used during prediction thresholding.
%      Default: 0.05
%
%  'relative_beta_SNR' :
%      Relative weighting between HbO and HbR misfit terms in the objective.
%      Larger values reduce the relative HbR penalty. Default: 1.0
%
%  'save_dir' :
%      Output directory for the saved absorption-estimation result MAT file.
%      Default: ''
%
%  'beta_od_file' :
%      Path to a ΔOD beta MAT file used to reconstruct Δμa if muA_file is
%      missing or empty. Default: ''
%
%  'save_recon_dir' :
%      Directory in which to save reconstructed Δμa if reconstruction is
%      triggered from opts.beta_od_file. Default: ''
%
%  'tik_alpha' :
%      Tikhonov regularization parameter used in Δμa reconstruction fallback.
%      Default: 0.01
%
%  'do_plots' :
%      If true, generate channel-level plots of estimated and measured HbO/HbR.
%      Default: true
%
%  'plot_clim' :
%      Color-axis limit used for optional channel plots. Default: 0.25
%
%  'plot_tag' :
%      String tag appended to output plot titles. Default: 'AbsorptionEstimation'
%
%  'plot_montage' :
%      Label passed through to the plotting function. Default: 'Generic'
%
%  'plot_thresh_meas' :
%      If true, also plot thresholded measured HbO/HbR maps. Default: true
%
%  Returns
%  -------
%  out : struct with fields
%
%    out.lambda_opt :
%        Selected regularization parameter from the L-curve.
%
%    out.lambda_grid :
%        Candidate lambda values evaluated during L-curve selection.
%
%    out.err_fit_vals :
%        Data-fit values across the lambda grid.
%
%    out.err_pen_vals :
%        Penalty-term values across the lambda grid.
%
%    out.x_opt :
%        Optimized hemoglobin parameter vector returned by the final fit.
%
%    out.dmu_a_opt :
%        Estimated spectral absorption magnitude vector:
%          [Δμa760; Δμa850]
%
%    out.errors :
%        Final objective components:
%          [total_error, HbO_error, HbR_error, ratio_penalty]
%
%    out.HbO_pred / out.HbR_pred :
%        Thresholded forward-predicted channel HbO/HbR values from the subject-
%        level Δμa reconstructions.
%
%    out.HbO_est / out.HbR_est :
%        Final fitted channel-level HbO/HbR estimates after threshold masking.
%
%    out.masks :
%        Struct containing channel-threshold masks:
%          .thres_HbO
%          .thres_HbR
%
%    out.SJ :
%        Spectral/chromophore-space Jacobian struct returned by
%        fnirspower.fwdcomp.buildPathAdjJacobian.
%
%    out.mean_mu_a760 / out.mean_mu_a850 :
%        Mean subject-level Δμa maps used as the spatial basis for fitting.
%
%    out.mu_a_est :
%        Struct containing fitted spectral absorption maps:
%          .mu_a760
%          .mu_a850
%
%  Example
%  -------
%  out = fnirspower.pipeline.run_absorption_estimation( ...
%      muA_file, nirsmodel_path, mesh_mat, beta_hb_file, layout_file, ...
%      'beta_od_file', beta_od_file, ...
%      'save_recon_dir', recon_dir, ...
%      'save_dir', out_dir, ...
%      'alpha', 0.05, ...
%      'relative_beta_SNR', 1.0, ...
%      'do_plots', true);
%
%  Example with Δμa reconstruction fallback
%  ----------------------------------------
%  out = fnirspower.pipeline.run_absorption_estimation( ...
%      '', nirsmodel_path, mesh_mat, beta_hb_file, layout_file, ...
%      'beta_od_file', beta_od_file, ...
%      'save_recon_dir', recon_dir, ...
%      'tik_alpha', 0.01);
%
%  Notes
%  -----
%  - If muA_file is missing or empty, opts.beta_od_file must be provided.
%  - The measured beta file must be compatible with
%    fnirspower.io.resolve_measured_betas.
%  - The plotting section auto-selects the first suitable variable in
%    layout_file; it does not assume a layout-specific variable name.
%  - This function performs absorption magnitude estimation, not instrument
%    calibration in the hardware sense.
%
%  See also
%  --------
%  fnirspower.recon.reconstruct_mu_a
%  fnirspower.fwdcomp.buildPathAdjJacobian
%  fnirspower.absmag.threshold_predictions
%  fnirspower.absmag.lcurve_lambda
%  fnirspower.absmag.muA_ls_obj_HbT_rel

arguments
    muA_file (1,:) char
    nirsmodel_path (1,:) char
    mesh_mat (1,:) char
    beta_file (1,:) char
    layout_file (1,:) char = ''
    opts.epsilon_mm_uM double = 2.303e-6 * [58.6,154.8; 105.8,69.1]
    opts.alpha (1,1) double = 0.05
    opts.relative_beta_SNR (1,1) double = 1.0
    opts.save_dir (1,:) char = ''
    opts.beta_od_file (1,:) char = ''
    opts.save_recon_dir (1,:) char = ''
    opts.tik_alpha (1,1) double = 0.01
    opts.do_plots (1,1) logical = true
    opts.plot_clim (1,1) double = 0.25
    opts.plot_tag (1,:) char = 'AbsorptionEstimation'
    opts.plot_montage (1,:) char = 'Generic'
    opts.plot_thresh_meas (1,1) logical = true
end

%% 0) Ensure Δμa exists (or reconstruct from ΔOD betas)
if ~(~isempty(muA_file) && exist(muA_file,'file')==2)
    if isempty(opts.beta_od_file) || exist(opts.beta_od_file,'file')~=2
        error(['No μa file found and opts.beta_od_file not provided. ', ...
            'Supply ΔOD betas (all_beta_a760/all_beta_a850) or pass an existing μa MAT.']);
    end
    
    muA_file = fnirspower.recon.reconstruct_mu_a( ...
        opts.beta_od_file, nirsmodel_path, mesh_mat, ...
        'save_dir', opts.save_recon_dir, ...
        'tik_alpha', opts.tik_alpha);
end

%% 1) Load Δμa data
S = load(muA_file, 'all_mu_a760','all_mu_a850','brain_nodes_idx');
if ~isfield(S,'all_mu_a760') || ~isfield(S,'all_mu_a850')
    error('muA_file must contain all_mu_a760 and all_mu_a850 (nSubj×nBrain).');
end

if ~isfield(S,'brain_nodes_idx') || isempty(S.brain_nodes_idx)
    M = load(mesh_mat,'genmesh');
    e = M.genmesh.ele;
    brain_tet = e(e(:,5)==2,:);
    S.brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);
end

mean_mu_a760 = mean(S.all_mu_a760, 1, 'omitnan')';
mean_mu_a850 = mean(S.all_mu_a850, 1, 'omitnan')';

%% 2) Spectral Jacobian (chromophore-space forward model)
nirs_model = load(nirsmodel_path);
SJ = fnirspower.fwdcomp.buildPathAdjJacobian(nirs_model, opts.epsilon_mm_uM, S.brain_nodes_idx);
% Expect SJ.J_chrom_brain sized: (2*nCh) × (2*nBrain)

%% 3) Load measured GLM betas
B = load(beta_file);
[meas_HbO, meas_HbR] = fnirspower.io.resolve_measured_betas(B);  % nCh×1 each
nCh = numel(meas_HbO);

%% 4) Forward-predict channel Hb per subject from Δμa
nSubj = size(S.all_mu_a760, 1);
all_Hb = NaN(nSubj, 2*nCh);
for i = 1:nSubj
    mu_vec = [S.all_mu_a760(i,:)'; S.all_mu_a850(i,:)'];  % (2*nBrain)×1
    all_Hb(i,:) = (SJ.J_chrom_brain * mu_vec).';
end
all_HbO_pred = all_Hb(:,1:nCh);
all_HbR_pred = all_Hb(:,nCh+1:end);

%% 5) Threshold predictions
[HbO_pred_thr, HbR_pred_thr, thres_HbO, thres_HbR] = ...
    fnirspower.absmag.threshold_predictions( ...
    meas_HbO, meas_HbR, all_HbO_pred, all_HbR_pred, opts.alpha);

meas_HbO_thr = meas_HbO;
meas_HbR_thr = meas_HbR;
meas_HbO_thr(thres_HbO) = 0;
meas_HbR_thr(thres_HbR) = 0;

%% 6) L-curve lambda selection (SNR-weighted fit term)
wR = 1 / max(eps, opts.relative_beta_SNR);

x0 = [5.70; -1.91];
lb = [0; -Inf];
ub = [Inf; 0];

obj_builder = @(lambda) deal( ...
    @(x) local_obj_scalar(x, lambda), ...
    @(x) local_eval_components(x, lambda));

[lambda_opt, lambda_grid, err_fit_vals, err_pen_vals] = ...
    fnirspower.absmag.lcurve_lambda(obj_builder, []);

%% 7) Final absorption estimation at lambda_opt
obj_final = @(x) local_obj_scalar(x, lambda_opt);
x_opt = fmincon(obj_final, x0, [], [], [], [], lb, ub, [], ...
    optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-8));

[err_tot, err_O, err_R, err_rat, dHb_opt, dmu_a_opt] = ...
    fnirspower.absmag.muA_ls_obj_HbT_rel( ...
    x_opt, mean_mu_a760, mean_mu_a850, thres_HbO, thres_HbR, ...
    opts.epsilon_mm_uM, SJ.J_chrom_brain, meas_HbO, meas_HbR, ...
    lambda_opt, opts.relative_beta_SNR);

%% 8) Channel-level prediction using fitted Δμa
m760 = mean_mu_a760 / max(abs(mean_mu_a760));
m850 = mean_mu_a850 / max(abs(mean_mu_a850));
vec_norm = [m760 * dmu_a_opt(1); m850 * dmu_a_opt(2)];
delta_Hb = SJ.J_chrom_brain * vec_norm;

HbO_est = delta_Hb(1:nCh);
HbR_est = delta_Hb(nCh+1:end);

HbO_est(thres_HbO) = 0;
HbR_est(thres_HbR) = 0;

%% 9) Pack outputs
out = struct();
out.lambda_opt   = lambda_opt;
out.lambda_grid  = lambda_grid;
out.err_fit_vals = err_fit_vals;
out.err_pen_vals = err_pen_vals;

out.x_opt     = dHb_opt;      % [ΔHbO; ΔHbR]
out.dmu_a_opt = dmu_a_opt;    % [Δμa760; Δμa850]
out.errors    = [err_tot, err_O, err_R, err_rat];

out.HbO_pred  = HbO_pred_thr;
out.HbR_pred  = HbR_pred_thr;
out.HbO_est   = HbO_est;
out.HbR_est   = HbR_est;

out.masks = struct('thres_HbO', thres_HbO, 'thres_HbR', thres_HbR);
out.SJ    = SJ;

out.mean_mu_a760 = mean_mu_a760;
out.mean_mu_a850 = mean_mu_a850;
out.mu_a_est.mu_a760 = mean_mu_a760 * dmu_a_opt(1);
out.mu_a_est.mu_a850 = mean_mu_a850 * dmu_a_opt(2);

%% 10) Optional plotting
if opts.do_plots && ~isempty(layout_file)
    L = load(layout_file);
    
    fn = fieldnames(L);
    layout2use = [];
    for k = 1:numel(fn)
        if isstruct(L.(fn{k})) || isnumeric(L.(fn{k}))
            layout2use = L.(fn{k});
            break;
        end
    end
    
    if ~isempty(layout2use) && exist('plot_channel_level_hb','file')==2
        date_str = datestr(datetime('now'),'yyyy-mm-dd');
        clim = opts.plot_clim;
        montage_type = opts.plot_montage;
        
        fnirspower.helpers.plot_channel_level_hb(layout2use, [HbO_est, HbR_est]', ...
            clim, date_str, [opts.plot_tag '_Estimated'], montage_type, 1);
        
        fnirspower.helpers.plot_channel_level_hb(layout2use, [meas_HbO', meas_HbR']', ...
            clim, date_str, [opts.plot_tag '_Measured'], montage_type, 1);
        
        if opts.plot_thresh_meas
            fnirspower.helpers.plot_channel_level_hb(layout2use, [meas_HbO_thr', meas_HbR_thr']', ...
                clim, date_str, [opts.plot_tag '_Measured_Thr'], montage_type, 1);
        end
    else
        warning('Plotting skipped: layout could not be resolved or plot_channel_level_hb not on path.');
    end
end

%% 11) Optional save
if ~isempty(opts.save_dir)
    if ~exist(opts.save_dir,'dir')
        mkdir(opts.save_dir);
    end
    date_str = datestr(datetime('now'),'yyyy-mm-dd');
    out_file = fullfile(opts.save_dir, [date_str '_Absorption_Results.mat']);
    save(out_file, 'out', '-v7.3');
    fprintf('Saved absorption estimation results to: %s\n', out_file);
end

%% ---- nested helpers ----
    function err = local_obj_scalar(x, lambda)
        [err, ~, ~, ~] = local_eval_components(x, lambda);
    end

    function [err, err_O, err_R, err_rat] = local_eval_components(x, lambda)
        [~, err_O, err_R_raw, err_rat] = ...
            fnirspower.absmag.muA_ls_obj_HbT_rel( ...
            x, mean_mu_a760, mean_mu_a850, ...
            thres_HbO, thres_HbR, ...
            opts.epsilon_mm_uM, ...
            SJ.J_chrom_brain, ...
            meas_HbO, meas_HbR, ...
            lambda, opts.relative_beta_SNR);
        
        % Return the weighted HbR error so the L-curve fit term matches the
        % data-fit term used by the optimization.
        err_R = wR * err_R_raw;
        
        err = err_O + err_R + lambda * err_rat;
        
        if ~isfinite(err)
            err = Inf;
        end
    end

end