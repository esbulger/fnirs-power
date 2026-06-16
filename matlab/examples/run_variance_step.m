%% RUN_VARIANCE_STEP
% Example driver for within-subject beta-variance estimation using
% fnirspower.pipeline.run_variance.
%
% This script:
%   1) sets up fnirspower and third-party paths
%   2) defines mesh / forward-model / layout inputs
%   3) builds an explicit SNIRF file list for the selected subjects
%   4) specifies block-resampling settings
%   5) runs the variance pipeline
%
% Notes
% -----
% - SNIRF files are resolved explicitly from workspace/rawdata so the run
%   does not depend on a fixed internal package path structure.
% - The most recent SNIRF file is selected if multiple files are found for
%   a subject.
% - This script currently uses a small subject subset and reduced iteration
%   settings for testing.

import fnirspower.*

%% -------------------- Setup paths and core inputs --------------------
P = fnirspower.setup_paths();

mesh_path    = fullfile(P.icbm_mesh_dir,'ICBM_mesh_5layer.mat');
nirsmodel_m  = fullfile(P.forward,'nirs_mesh_ICBM_5layer_VisualDense_nirsmodel.mat');
layout_file  = fullfile(P.layouts,'layout_VisualDense.mat');

%% -------------------- Subject selection --------------------
% Full subject list:
% subjects = int16([101 102 103 105 106 107 108 110 111 112 113 114 116 117 118 119]);

% Small subset for test runs:
subjects = int16([101 102]);

%% -------------------- Output location --------------------
save_dir = fullfile(P.derivatives,'variance');
if exist(save_dir,'dir') ~= 7
    mkdir(save_dir);
end

%% -------------------- Resolve SNIRF files --------------------
RAW_ROOT = fullfile(P.workspace,'rawdata');

% Build explicit SNIRF file list in the same order as subjects.
snirf_files = cell(numel(subjects),1);
for i = 1:numel(subjects)
    sid = subjects(i);
    snirf_dir = fullfile(RAW_ROOT, sprintf('Subject%d', sid), 'leftright');
    d = dir(fullfile(snirf_dir, '*.snirf'));
    assert(~isempty(d), 'No SNIRF found for subject %d under %s', sid, snirf_dir);

    % Choose most recent if multiple files are present.
    [~,ix] = sort([d.datenum], 'descend');
    snirf_files{i} = fullfile(d(ix(1)).folder, d(ix(1)).name);
end

paths = struct();
paths.snirf_files = snirf_files;

%% -------------------- Variance-analysis settings --------------------
baseline_sec = 2.5;
block_sec    = 15;
iter_blocks  = [4 8 12 16];
n_iter       = 10;

% Use the same baseline duration when concatenating resampled blocks.
concat_baseline_sec = baseline_sec;

% Plotting controls
do_plots     = true;

% Date tag used for saved outputs
date_str = datestr(datetime('now'),'yyyy-mm-dd');

%% -------------------- Run summary --------------------
fprintf('[run_variance_step]\n');
fprintf('  n SNIRFs  : %d\n', numel(paths.snirf_files));
fprintf('  model     : %s\n', nirsmodel_m);
fprintf('  mesh      : %s\n', mesh_path);
fprintf('  save_dir  : %s\n', save_dir);

%% -------------------- Run pipeline --------------------
V = fnirspower.pipeline.run_variance( ...
    double(subjects), paths, mesh_path, nirsmodel_m, layout_file, ...
    'baseline_sec', baseline_sec, ...
    'block_sec', block_sec, ...
    'iter_blocks', iter_blocks, ...
    'n_iter', n_iter, ...
    'concat_baseline_sec', concat_baseline_sec, ...
    'do_plots', do_plots, ...
    'layout_var', '', ...
    'save_dir', save_dir, ...
    'date_str', date_str);

fprintf('[run_variance_step] Done.\n');