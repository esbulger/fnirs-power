%% RUN_GLM_STEP (group orchestrator version)
% Uses pipeline.run_group_glm (which calls pipeline.run_subject_glm) to:
%  - preprocess + GLM per subject (HbO/HbR and ΔOD optional)
%  - collect per-subject arrays
%  - save the same files your previous script produced
%
%  Outputs to SAVE_DIR:
%    <DATE>_GLM_Betas_Hb.mat            (subjects, all_subj_beta_hbo/hbr, mean_hbo/hbr)
%    <DATE>_GLM_Betas_Hb_Masked.mat     (Visual_Checkerboard_GLM_Mean_Subj_Beta_Channel_Masked)
%    <DATE>_GLM_Betas_OD.mat            (subjects, all_beta_a760/all_beta_a850)
%
% Plots (optional) replicate your previous behavior.

%% --- Resolve project-relative paths ---
import fnirspower.*
% P = fnirspower.setup_paths();

% Inputs (relative)
raw_root    = fullfile(P.workspace,'rawdata');                 % SubjectXXX/leftright/*.snirf
mesh_path   = fullfile(P.icbm_mesh_dir,'ICBM_mesh_5layer.mat');% contains genmesh.node/ele
nirs_path = fullfile(P.forward,'nirs_mesh_ICBM_5layer_VisualDense_nirsmodel.mat');
layout_path  = fullfile(P.layouts,'layout_VisualDense.mat');
montage_type = 'Visual_Dense';

% Outputs
save_dir = fullfile(P.derivatives,'glm');  % e.g., workspace/derivatives/glm
if ~exist(save_dir,'dir'), mkdir(save_dir); end

% Subject list
subjects = [101 102 103 105 106 107 108 110 111 112 113 114 116 117 118 119];
% subjects = [101 102];

% Flags
DO_PLOTS = true;

%% -------------------- Run group GLM (delegated) --------------------
% Build explicit SNIRF file list in the same order as subjects
snirf_files = cell(numel(subjects),1);
for i = 1:numel(subjects)
    sid = subjects(i);
    snirf_dir = fullfile(raw_root, sprintf('Subject%d', sid), 'leftright');
    d = dir(fullfile(snirf_dir, '*.snirf'));
    assert(~isempty(d), 'No SNIRF found for subject %d under %s', sid, snirf_dir);
    
    % choose most recent if multiple
    [~,ix] = sort([d.datenum], 'descend');
    snirf_files{i} = fullfile(d(ix(1)).folder, d(ix(1)).name);
end

snirf_paths = struct();
snirf_paths.snirf_files = snirf_files;

% Run group level glm
G = fnirspower.pipeline.run_group_glm( ...
    subjects, snirf_paths, mesh_path, nirs_path, layout_path, ...
    'hrf_seconds', 32, ...
    'alpha', 0.01, ...
    'n_pca', 1, ...
    'r2max', 2.5, ...
    'Pmax', []);

% Convenience names for saving (match your old filenames/vars)
all_subj_beta_hbo = G.beta.hbo;
all_subj_beta_hbr = G.beta.hbr;
all_beta_a760     = G.beta.od760;   % ΔOD betas @760
all_beta_a850     = G.beta.od850;   % ΔOD betas @850

% Means (same as before)
mean_hbo = mean(all_subj_beta_hbo, 1, 'omitnan');
mean_hbr = mean(all_subj_beta_hbr, 1, 'omitnan');

% Group p-masking
[~, p_hbo] = ttest(all_subj_beta_hbo, 0);
[~, p_hbr] = ttest(all_subj_beta_hbr, 0);
mask_hbo = p_hbo > 0.05;
mask_hbr = p_hbr > 0.05;
mean_hbo_m = mean_hbo;
mean_hbr_m = mean_hbr;
mean_hbo_m(mask_hbo) = 0;
mean_hbr_m(mask_hbr) = 0;

%% -------------------- Save --------------------
date_str = char(datetime('now', 'Format', 'yyyy-MM-dd'));

% Hb betas: per-subject, per-channel values and group means.
GLM_Hb = struct();
GLM_Hb.subjects = subjects;
GLM_Hb.all_subj_beta_hbo = all_subj_beta_hbo;
GLM_Hb.all_subj_beta_hbr = all_subj_beta_hbr;
GLM_Hb.mean_hbo = mean_hbo;
GLM_Hb.mean_hbr = mean_hbr;

hb_file = fullfile( ...
    save_dir, ...
    sprintf('%s_GLM_Betas_Hb.mat', date_str));

save(hb_file, 'GLM_Hb', '-v7.3');

% Masked group-mean Hb betas.
GLM_Hb_Masked = struct();
GLM_Hb_Masked.mean_hbo = mean_hbo_m;
GLM_Hb_Masked.mean_hbr = mean_hbr_m;

hb_masked_file = fullfile( ...
    save_dir, ...
    sprintf('%s_GLM_Betas_Hb_Masked.mat', date_str));

save(hb_masked_file, 'GLM_Hb_Masked', '-v7.3');

% Delta-OD betas used for absorption reconstruction.
GLM_OD = struct();
GLM_OD.subjects = subjects;
GLM_OD.all_beta_a760 = all_beta_a760;
GLM_OD.all_beta_a850 = all_beta_a850;

od_file = fullfile( ...
    save_dir, ...
    sprintf('%s_GLM_Betas_OD.mat', date_str));

save(od_file, 'GLM_OD', '-v7.3');

%% ---- Also save μa recon so dmua can load directly ----
% G.mua.lambda760 / lambda850 are nSubj × nBrain (brain-only nodes)
all_mu_a760     = G.mua.lambda760;
all_mu_a850     = G.mua.lambda850;

% Recover brain_nodes_idx from the mesh passed into run_group_glm
M = load(mesh_path,'genmesh');
e = M.genmesh.ele;
brain_tet = e(e(:,5)==2,:);
brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);

muA_fn = fullfile(fullfile(P.derivatives,'recon'), sprintf('%s_mu_a_recon.mat', date_str));
save(muA_fn, 'all_mu_a760','all_mu_a850','brain_nodes_idx','-v7.3');
fprintf('Saved μa recon for absorption estimation: %s\n', muA_fn);

%% -------------------- Plots (optional) --------------------
if DO_PLOTS
    L = load(layout_path);
    fL = fieldnames(L);
    assert(~isempty(fL), 'No variables found in layout file: %s', layout_path);
    layout2use = L.(fL{1});
    
    % Group average
    meanBetaGLM = [mean_hbo; mean_hbr];
    clim1 = max(max(abs(meanBetaGLM)));
    plot_channel_level_hb(layout2use, meanBetaGLM, ...
        clim1, date_str, 'GLM Block Avg Beta', montage_type, 1);
    
    % Uncorrected p-masked
    clim2 = max(max(abs([mean_hbo_m; mean_hbr_m])));
    plot_channel_level_hb(layout2use, [mean_hbo_m; mean_hbr_m], ...
        clim2, date_str, 'GLM Block Avg Beta (p-masked)', montage_type, 1);
    
    
    %% ---------- Group-level source reconstruction plots ----------
    Mmesh = load(mesh_path);
    if isfield(Mmesh, 'genmesh')
        p = Mmesh.genmesh.node(:,1:3);
        e = Mmesh.genmesh.ele;
    else
        p = Mmesh.p;
        e = Mmesh.e;
    end
    
    % Mean reconstructions across subjects
    mean_recon_hbo = mean(G.recon.hbo, 1, 'omitnan');
    mean_recon_hbr = mean(G.recon.hbr, 1, 'omitnan');
    
    % HbO reconstruction
    figure('Position',[600 200 1200 800]); hold on;
    iso2mesh_plotmesh([p mean_recon_hbo'], e(e(:,5)==2,:), 'EdgeAlpha', 0.2, 'FaceAlpha', 1);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(['Mean HbO Reconstruction ' date_str]);
    colormap(redblue);
    colorbar;
    view(0,15);
    clim_hbo = max(abs(mean_recon_hbo'));
    if clim_hbo > 0
        caxis([-clim_hbo, clim_hbo]);
    end
    
    % HbR reconstruction
    figure('Position',[600 200 1200 800]); hold on;
    iso2mesh_plotmesh([p mean_recon_hbr'], e(e(:,5)==2,:), 'EdgeAlpha', 0.2, 'FaceAlpha', 1);
    xlabel('x'); ylabel('y'); zlabel('z');
    title(['Mean HbR Reconstruction ' date_str]);
    colormap(redblue);
    colorbar;
    view(0,15);
    clim_hbr = max(abs(mean_recon_hbr));
    if clim_hbr > 0
        caxis([-clim_hbr, clim_hbr]);
    end
end