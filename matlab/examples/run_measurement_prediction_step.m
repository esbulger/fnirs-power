%% RUN_MEASUREMENT_PREDICTION_STEP
%  Driver: resolve project paths, pick ROI(s), run measurement prediction,
%  optionally plot.

import fnirspower.*

P = fnirspower.setup_paths();

%% Load head model, forward model, dmu_a information
montageName = 'dense';

% Files
mesh_mat       = fullfile(P.icbm_mesh_dir,'ICBM_mesh_5layer.mat');
nirsmodel_file = fullfile(P.forward,['nirs_mesh_ICBM_5layer_', montageName, '_nirsmodel.mat']);
layout_file = fullfile(P.layouts, 'layout_dense.mat');

% plotting options
plotGreyMatter = 0;
plotROIOnMesh = 0;

%% Estimated / assumed spectral amplitude
% ----------------------------------------
% If prediction_mode = 'estimated':
%   - set dmua_file explicitly based on absorption estimation step
%
%   OR
%
% If prediction_mode = 'assumed':
%   - dmua_file is ignored
%   - dmu_a_assumed = [dmu760; dmu850] is used directly

prediction_mode = 'estimated';   % 'estimated' or 'assumed'

dmua_file = fullfile(P.derivatives, 'absmag', 'Default_Absorption_Results.mat');
dmu_a_assumed = [1e-6; 1e-4];    % used only if prediction_mode = 'assumed';

%% Optional grey-matter mesh plot

if plotGreyMatter
    mesh_model = load(mesh_mat);

    if isfield(mesh_model,'genmesh')
        p = mesh_model.genmesh.node(:,1:3);
        e = mesh_model.genmesh.ele;
    else
        p = mesh_model.p;
        e = mesh_model.e;
    end

    figure('Position',[600 200 1200 800]); hold on;
    iso2mesh_plotmesh(p, e(e(:,5)==1,:), 'EdgeAlpha', 0.2, 'FaceColor', [1, 1, 1], 'FaceAlpha', 1);
    iso2mesh_plotmesh(p, e(e(:,5)==2,:), 'EdgeAlpha', 0.2, 'FaceColor', [1, 1, 1], 'FaceAlpha', 1);
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Brain mesh');
    view([90, 0])
    axis equal;
    axis off;
end


%% ROI definition
roi_def = struct();
roi_def.centers_mm = [70, -10, 9; ...
                      70, -10, 11];   % 1 ROI (make this [nROI x 3] for multiple ROIs)
roi_def.radius_mm  = [7.5, 7.5];
roi_def.sigma_mm   = [5, 5];
roi_def.amp_std        = [0.3, 0.3];
roi_def.names      = {'ROI', 'ROI2'};  % optional, dimension must match above

%% Optional ROI plot on grey-matter mesh
if plotROIOnMesh
    mesh_model = load(mesh_mat);

    if isfield(mesh_model,'genmesh')
        p = mesh_model.genmesh.node(:,1:3);
        e = mesh_model.genmesh.ele;
    else
        p = mesh_model.p;
        e = mesh_model.e;
    end

    gm_ele = e(e(:,5)==2,:);
    gm_nodes_idx = unique(gm_ele(:,1:4));
    p_gm = p(gm_nodes_idx,:);

    roi_center = roi_def.centers_mm(1,:);
    roi_radius = roi_def.radius_mm(1);

    d = sqrt(sum((p_gm - roi_center).^2, 2));
    roi_mask = d <= roi_radius;

    roi_overlay = zeros(size(p,1),1);
    roi_overlay(gm_nodes_idx(roi_mask)) = 1;

    p_plot = [p roi_overlay];

    figure('Position',[600 200 1200 800]); hold on;
    iso2mesh_plotmesh(p, e(e(:,5)==1,:), 'EdgeAlpha', 0.2, 'FaceColor', [1, 1, 1], 'FaceAlpha', 1);
    iso2mesh_plotmesh(p_plot, gm_ele, 'EdgeAlpha', 0.2, 'FaceAlpha', 1);

    colormap(redblue);
    caxis([-1 1]);

    xlabel('x'); ylabel('y'); zlabel('z');
    view([90, 0])
    axis equal;
    axis off;
end

%% Common options
nSubjects        = 50;
hbRatio          = 3;
plotMesh         = true;
plotROISeed      = true;
plotAvgProb      = true;
plotSubjects     = 0;
computeDistances = true;
plotChannels     = true;

%% Run prediction function
switch lower(prediction_mode)
    case 'estimated'
        M = fnirspower.pipeline.run_measurement_prediction( ...
            mesh_mat, nirsmodel_file, roi_def, dmua_file, ...
            'nSubjects', nSubjects, ...
            'hbRatio', hbRatio, ...
            'dmuaScale', true, ...
            'dmuaScaleDiv', 3, ...  % 3% for visual checkerboard stimuli
            'plotMesh', plotMesh, ...
            'plotROISeed', plotROISeed, ...
            'plotAvgProb', plotAvgProb, ...
            'plotSubjects', plotSubjects, ...
            'computeDistances', computeDistances, ...
            'plotChannels', plotChannels, ...
            'layout_file', layout_file, ...
            'montageName', montageName);

    case 'assumed'
        M = fnirspower.pipeline.run_measurement_prediction( ...
            mesh_mat, nirsmodel_file, roi_def, '', ...
            'dmu_a', dmu_a_assumed, ...
            'dmuaScale', false, ...
            'nSubjects', nSubjects, ...
            'hbRatio', hbRatio, ...
            'plotMesh', plotMesh, ...
            'plotROISeed', plotROISeed, ...
            'plotAvgProb', plotAvgProb, ...
            'plotSubjects', plotSubjects, ...
            'computeDistances', computeDistances, ...
            'plotChannels', plotChannels, ...
            'layout_file', layout_file, ...
            'montageName', montageName);

    otherwise
        error('prediction_mode must be ''estimated'' or ''assumed''.');
end

%% Save
outdir = fullfile(P.derivatives,'predict');
if exist(outdir,'dir')~=7, mkdir(outdir); end
date_str = datestr(datetime('now'),'yyyy-mm-dd');
outfile = fullfile(outdir, sprintf('%s_measurement_prediction.mat', date_str));
save(outfile, 'M', '-v7.3');
fprintf('Saved: %s\n', outfile);
