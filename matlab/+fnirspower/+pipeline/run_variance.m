function V = run_variance(subject_ids, paths, mesh_path, nirsmodel_path, layout_mat, opts)
%RUN_VARIANCE Estimate within-subject GLM beta variability using block-resampling.
%
%  V = fnirspower.pipeline.run_variance( ...
%      subject_ids, paths, mesh_path, nirsmodel_path, layout_mat, ...)
%
%  Description
%  -----------
%  This function estimates within-subject variability in channel-level GLM
%  beta values by repeatedly resampling task blocks and refitting subject-
%  level GLMs. It:
%
%    1) Loads the shared mesh and forward-model information used for analysis.
%    2) Resolves one SNIRF file per subject from paths.snirf_files.
%    3) Runs subject-level block-resampling analysis through
%       fnirspower.pipeline.run_variance_subject.
%    4) Aggregates HbO/HbR beta estimates across subjects, iterations, and
%       block-count conditions.
%    5) Computes summary variability curves as a function of block count.
%    6) Fits a simple A/sqrt(B) model to the variability curve.
%    7) Optionally plots the summary curve.
%    8) Optionally saves the full result structure to disk.
%
%  Required inputs
%  ---------------
%  subject_ids :
%      Numeric vector of subject IDs. The order of subject_ids must match
%      the order of paths.snirf_files.
%
%  paths :
%      Struct containing:
%
%        paths.snirf_files :
%            Cell array or string array of full SNIRF file paths, one per
%            subject, in the same order as subject_ids.
%
%  mesh_path :
%      Path to the tetrahedral head mesh MAT file. Used to recover brain-node
%      indices and populate the model struct passed to the subject-level
%      variance routine.
%
%  nirsmodel_path :
%      Path to the forward-model MAT file used during subject-level analysis.
%
%  layout_mat :
%      Path to a layout MAT file for optional topographic plotting.
%      Pass '' to skip layout loading.
%
%  Name-value options
%  ------------------
%  'baseline_sec' :
%      Baseline duration in seconds used during block extraction and
%      concatenation. Default: 2.5
%
%  'block_sec' :
%      Duration of each task block in seconds. Default: 15
%
%  'iter_blocks' :
%      Vector of block counts to evaluate during resampling.
%      Default: [4 8 12 16 20 24]
%
%  'n_iter' :
%      Number of resampling iterations per block-count condition.
%      Default: 100
%
%  'concat_baseline_sec' :
%      Baseline duration used when concatenating resampled blocks.
%      If left unspecified, this defaults to baseline_sec.
%      Default: baseline_sec
%
%  'do_plots' :
%      If true, plot the summary variability curve and fitted A/sqrt(B)
%      relationship. Default: false
%
%
%  'layout_var' :
%      Optional variable name to load from layout_mat. If empty, the first
%      variable in the layout MAT file is used. Default: ''
%
%  'save_dir' :
%      Directory to which the full result MAT file should be saved.
%      Default: ''
%
%  'date_str' :
%      Date tag used in optional saved output filenames and topographic plot
%      labels. Default: datestr(datetime('now'),'yyyy-mm-dd')
%
%  Returns
%  -------
%  V : struct with fields
%
%    V.subjects :
%        Subject IDs included in the run.
%
%    V.iter_blocks :
%        Vector of block counts evaluated.
%
%    V.n_iter :
%        Number of bootstrap/resampling iterations per condition.
%
%    V.beta_hbo :
%        HbO beta estimates with size:
%          [nSubj x nIter x nBlockConds x nCh]
%
%    V.beta_hbr :
%        HbR beta estimates with size:
%          [nSubj x nIter x nBlockConds x nCh]
%
%    V.summary :
%        Struct containing summary variability metrics:
%
%          .std_iter_hbo :
%              Standard deviation across iterations for each subject,
%              block-count condition, and channel.
%
%          .curve_mean :
%              Mean variability curve collapsed across subjects/channels.
%
%          .curve_std :
%              Standard deviation of the variability curve across
%              subject-channel combinations.
%
%          .fit_A :
%              Scalar coefficient of the fitted A/sqrt(B) model.
%
%          .fit_y :
%              Fitted curve values at iter_blocks.
%
%          .fit_R2 :
%              Coefficient of determination for the fitted A/sqrt(B) model.
%
%
%    V.layout :
%        Raw loaded layout struct, if layout_mat was provided and found.
%
%  Example
%  -------
%  V = fnirspower.pipeline.run_variance( ...
%      subject_ids, paths, mesh_path, nirsmodel_path, layout_mat, ...
%      'baseline_sec', 2.5, ...
%      'block_sec', 15, ...
%      'iter_blocks', [4 8 12 16 20 24], ...
%      'n_iter', 100, ...
%      'concat_baseline_sec', 2.5, ...
%      'do_plots', true, ...
%      'layout_var', '', ...
%      'save_dir', save_dir, ...
%      'date_str', datestr(datetime('now'),'yyyy-mm-dd'));
%
%  Notes
%  -----
%  - subject_ids and paths.snirf_files must be aligned in the same order.
%  - This function relies on fnirspower.pipeline.run_variance_subject for
%    subject-level preprocessing, block resampling, and GLM fitting.
%  - Topographic plotting requires:
%      1) a valid layout MAT file,
%      2) an appropriate layout variable inside that MAT file, and
%      3) fnirspower.helpers.plot_channel_level_hb on the MATLAB path.
%  - If one or more SNIRF files are missing, those subjects are skipped in
%    the subject loop; an error is only raised if no subjects produce output.
%
%  See also
%  --------
%  fnirspower.pipeline.run_variance_subject
%  fnirspower.variance.epoch_blocks
%  fnirspower.variance.concat_blocks
%  fnirspower.helpers.fnirspower.helpers.plot_channel_level_hb

arguments
    subject_ids (1,:) double
    paths struct
    mesh_path (1,:) char
    nirsmodel_path (1,:) char
    layout_mat (1,:) char = ''
    opts.baseline_sec (1,1) double = 2.5
    opts.block_sec (1,1) double = 15
    opts.iter_blocks (1,:) double = [4 8 12 16 20 24]
    opts.n_iter (1,1) double = 100
    opts.concat_baseline_sec (1,1) double = NaN
    opts.do_plots (1,1) logical = false
    opts.layout_var (1,:) char = ''
    opts.save_dir (1,:) char = ''
    opts.date_str (1,:) char = datestr(datetime('now'),'yyyy-mm-dd')
end

if isnan(opts.concat_baseline_sec)
    opts.concat_baseline_sec = opts.baseline_sec;
end

if ~isfield(paths, 'snirf_files') || isempty(paths.snirf_files)
    error('paths.snirf_files is required and must match subject_ids.');
end

snirf_files = paths.snirf_files;
if isstring(snirf_files)
    snirf_files = cellstr(snirf_files(:));
end
if ~iscell(snirf_files)
    error('paths.snirf_files must be a cell array or string array of full SNIRF paths.');
end
if numel(snirf_files) ~= numel(subject_ids)
    error('paths.snirf_files must have the same length/order as subject_ids.');
end

%% ---------------- Load mesh / model once ----------------
disp("Loading mesh...")
M = load(mesh_path);
if isfield(M,'genmesh')
    e = M.genmesh.ele;
else
    e = M.e;
end
brain_tet = e(e(:,5)==2,:);
brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);

NM = load(nirsmodel_path);
model = NM;
model.brain_nodes_idx = brain_nodes_idx;

if ~isfield(model,'epsilon_mm_uM')
    model.epsilon_mm_uM = 2.303e-6 * [58.6,154.8; 105.8,69.1];
end

if ~isempty(layout_mat) && exist(layout_mat,'file')==2
    L = load(layout_mat);
else
    L = [];
end

%% ---------------- Per-subject run (PARFOR-safe) ----------------
disp("Initializing parallel pool for subject-level looping...")
nSubj = numel(subject_ids);
iter_blocks = opts.iter_blocks(:).';
nB = numel(iter_blocks);

Ycell = cell(nSubj,1);

parfor i = 1:nSubj
    fprintf("Beginning variance estimation for subject %i", i)
    sid = subject_ids(i);
    snirf_path = char(snirf_files{i});

    if exist(snirf_path,'file') ~= 2
        warning('SNIRF file not found for subject %d: %s', sid, snirf_path);
        Ycell{i} = [];
        continue;
    end

    OUT = fnirspower.pipeline.run_variance_subject( ...
        snirf_path, model, opts.baseline_sec, opts.block_sec, ...
        iter_blocks, opts.n_iter, opts.concat_baseline_sec);

    Ycell{i} = OUT;
end

%% ---------------- Pack into arrays ----------------
first_ok = find(~cellfun(@isempty,Ycell), 1, 'first');
if isempty(first_ok)
    error('run_variance: no subjects produced output.');
end

nCh = size(Ycell{first_ok}.beta_hbo, 3);
beta_hbo_all = NaN(nSubj, opts.n_iter, nB, nCh);
beta_hbr_all = NaN(nSubj, opts.n_iter, nB, nCh);

for i = 1:nSubj
    if isempty(Ycell{i}), continue; end
    beta_hbo_all(i,:,:,:) = Ycell{i}.beta_hbo;
    beta_hbr_all(i,:,:,:) = Ycell{i}.beta_hbr;
end

%% ---------------- Summary metrics ----------------
std_iter_hbo = squeeze(std(beta_hbo_all, 0, 2, 'omitnan'));  % [nSubj x nB x nCh]
std_rs = reshape(permute(std_iter_hbo, [2 1 3]), [nB, nSubj*nCh]).';

mean_curve = mean(std_rs, 1, 'omitnan');
std_curve  = std(std_rs, 0, 1, 'omitnan');

fit_out = fnirspower.variance.fit_A_over_sqrt_blocks(iter_blocks, mean_curve);
A_fit = fit_out.A;
y_fit = fit_out.y_fit;
R2 = fit_out.R2;

summary = struct();
summary.std_iter_hbo = std_iter_hbo;
summary.curve_mean   = mean_curve;
summary.curve_std    = std_curve;
summary.fit_A        = A_fit;
summary.fit_y        = y_fit(:).';
summary.fit_R2       = R2;

%% ---------------- Fit summary printout ----------------
fprintf('\n[run_variance] Variance summary\n');
fprintf('  Subjects included      : %d\n', nSubj);
fprintf('  Iterations per block   : %d\n', opts.n_iter);
fprintf('  Block counts evaluated : %s\n', mat2str(iter_blocks));
fprintf('  Channels analyzed      : %d\n', nCh);
fprintf('  Subject-channel pairs  : %d\n', nSubj * nCh);
fprintf('  Mean std(beta_HbO)     : %s\n', mat2str(mean_curve, 4));
fprintf('  Std of std(beta_HbO)   : %s\n', mat2str(std_curve, 4));
fprintf(['  Fitted relationship    : std(beta_HbO) = A / sqrt(B), ' ...
         'with A = %.4f and R^2 = %.4f\n'], A_fit, R2);
fprintf(['  Interpretation         : the expected within-subject HbO beta ' ...
         'variability decreases approximately as 1/sqrt(B), where B is the number of blocks.\n\n']);

%% ---------------- Optional plot of fit ----------------
if opts.do_plots
    figure; hold on;
    errorbar(iter_blocks, mean_curve, std_curve, '-o');
    plot(iter_blocks, summary.fit_y, '--');
    xlabel('Block Count');
    ylabel('Bootstrapped Std (HbO Beta)');
    title(sprintf('Within-subject Beta variability: A/sqrt(B), A=%.3f, R^2=%.3f', A_fit, R2));
    grid on;
end

%% ---------------- Output struct ----------------
V = struct();
V.subjects    = subject_ids;
V.iter_blocks = iter_blocks;
V.n_iter      = opts.n_iter;
V.beta_hbo    = beta_hbo_all;
V.beta_hbr    = beta_hbr_all;
V.summary     = summary;
if ~isempty(L), V.layout = L; end

%% ---------------- Optional save ----------------
if ~isempty(opts.save_dir)
    if exist(opts.save_dir,'dir')~=7
        mkdir(opts.save_dir);
    end
    out_fn = fullfile(opts.save_dir, sprintf('%s_Variance_HbO.mat', opts.date_str));
    save(out_fn, 'V', '-v7.3');
    fprintf('Saved variance results: %s\n', out_fn);
end

end
