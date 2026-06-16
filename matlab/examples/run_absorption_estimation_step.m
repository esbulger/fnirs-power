%% RUN_ABSORPTION_ESTIMATION_STEP
% Driver script: ensure μa exists, then run absorption estimation.
% - Finds latest GLM outputs (Hb and ΔOD) under project paths.
% - If μa recon is missing, reconstructs it from ΔOD betas.
% - Calls fnirspower.pipeline.run_absorption_estimation with resolved files.

%% Resolve project paths
P = fnirspower.paths();
addpath(P.matlab);

%% User options
use_masked = true;
muA_file = '';
beta_od_file = '';
beta_hb_file = '';

% NIRFAST model + mesh + layout
nirsmodel_file = fullfile(P.forward, 'nirs_mesh_ICBM_5layer_VisualDense_nirsmodel.mat');
mesh_mat       = fullfile(P.icbm_mesh_dir, 'ICBM_mesh_5layer.mat');
layout_file    = fullfile(P.layouts, 'layout_VisualDense.mat');

%% Locate latest GLM Hb (new or legacy)
if isempty(beta_hb_file)
    hb_new = dir(fullfile(P.derivatives, 'glm', '*_GLM_Betas_Hb.mat'));
    if ~isempty(hb_new)
        [~,ix] = sort([hb_new.datenum], 'descend');
        beta_hb_file = fullfile(hb_new(ix(1)).folder, hb_new(ix(1)).name);
    else
        hb_legacy = dir(fullfile(P.derivatives, 'glm', '*_GLM_Betas_Hb_Masked.mat'));
        if isempty(hb_legacy)
            error('No GLM Hb files found under %s', fullfile(P.derivatives, 'glm'));
        end
        [~,ix] = sort([hb_legacy.datenum], 'descend');
        beta_hb_file = fullfile(hb_legacy(ix(1)).folder, hb_legacy(ix(1)).name);
    end
end

%% Prefer masked Hb means when requested
if use_masked
    mlist = dir(fullfile(P.derivatives, 'glm', '*_GLM_Betas_Hb_Masked.mat'));
    if ~isempty(mlist)
        [~,ix] = sort([mlist.datenum], 'descend');
        beta_hb_file = fullfile(mlist(ix(1)).folder, mlist(ix(1)).name);
    end
end

%% Locate latest ΔOD betas (needed if μa recon must be built)
if isempty(beta_od_file)
    odlist = dir(fullfile(P.derivatives, 'glm', '*_GLM_Betas_OD.mat'));
    if ~isempty(odlist)
        [~,ix] = sort([odlist.datenum], 'descend');
        beta_od_file = fullfile(odlist(ix(1)).folder, odlist(ix(1)).name);
    else
        beta_od_file = '';
    end
end

%% Locate existing μa recon if present
if isempty(muA_file)
    rlist = dir(fullfile(P.derivatives, 'recon', '*_mu_a_recon.mat'));
    if ~isempty(rlist)
        [~,ix] = sort([rlist.datenum], 'descend');
        muA_file = fullfile(rlist(ix(1)).folder, rlist(ix(1)).name);
    else
        muA_file = '';
    end
end

%% Absorption estimation options
epsilon_mm_uM = 2.303e-6 * [58.6,154.8; 105.8,69.1];
alpha = 0.05;
relative_beta_SNR = fnirspower.io.try_load_relative_snr(fullfile(P.derivatives, 'glm'));
if isnan(relative_beta_SNR)
    relative_beta_SNR = 1.0;
end

save_dir = fullfile(P.derivatives, 'absmag');
save_recon_dir = fullfile(P.derivatives, 'recon');
tik_alpha = 0.01;

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end
if exist(save_recon_dir, 'dir') ~= 7
    mkdir(save_recon_dir);
end

do_plots = true;
plot_clim = 0.25;
plot_tag = 'Visual_Checkerboard';
plot_montage = 'Visual_Dense';

fprintf(['[run_absorption_estimation_step]\n' ...
    '  Hb file   : %s\n' ...
    '  OD file   : %s\n' ...
    '  μa file   : %s\n' ...
    '  model     : %s\n' ...
    '  mesh      : %s\n'], ...
    beta_hb_file, ...
    beta_od_file, ...
    fnirspower.helpers.ternary(~isempty(muA_file), muA_file, '<build>'), ...
    nirsmodel_file, ...
    mesh_mat);

%% Execute absorption estimation
% If muA_file is empty, it will be reconstructed using beta_od_file.
out = fnirspower.pipeline.run_absorption_estimation( ...
    muA_file, nirsmodel_file, mesh_mat, beta_hb_file, layout_file, ...
    'epsilon_mm_uM', epsilon_mm_uM, ...
    'alpha', alpha, ...
    'relative_beta_SNR', relative_beta_SNR, ...
    'save_dir', save_dir, ...
    'beta_od_file', beta_od_file, ...
    'save_recon_dir', save_recon_dir, ...
    'tik_alpha', tik_alpha, ...
    'do_plots', do_plots, ...
    'plot_clim', plot_clim, ...
    'plot_tag', plot_tag, ...
    'plot_montage', plot_montage, ...
    'plot_thresh_meas', true);

fprintf('[run_absorption_estimation_step] Done.\n');