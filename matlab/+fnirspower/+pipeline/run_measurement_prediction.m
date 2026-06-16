function M = run_measurement_prediction(mesh_mat, nirsmodel_file, roi_def, dmua_file, opts)
%RUN_MEASUREMENT_PREDICTION Predict noiseless channel-level responses from simulated ROI distributions.
%
%  M = fnirspower.pipeline.run_measurement_prediction( ...
%      mesh_mat, nirsmodel_file, roi_def, dmua_file, ...)
%
%  Description
%  -----------
%  This function generates subject-level probabilistic ROI maps, converts
%  them into voxelwise spectral absorption-change maps, and
%  forward-projects those maps to the obtain channel level dHbO/dHbR.
%  It is intended for measurement prediction and ROI-based simulation.
%
%  The function:
%    1) Loads the tetrahedral mesh and extracts brain nodes.
%    2) Loads the forward model and derives channel midpoint locations.
%    3) Resolves a spectral absorption magnitude either from a saved
%       absorption-estimation result file or from a directly supplied
%       [dmu760; dmu850] vector.
%    4) Generates subject-level probabilistic ROI maps for one or more ROIs.
%    5) Optionally visualizes the mesh, ROI seed(s), average ROI maps, and
%       example subject maps using iso2mesh-style plotting.
%    6) Converts ROI probability maps into voxelwise HbO/HbR and Δμa maps.
%    7) Forward-projects those maps to the channel level using the spectral
%       Jacobian.
%    8) Optionally visualizes channel-level predicted HbO/HbR maps.
%    9) Optionally computes distances between channel midpoints and ROI
%       centroids.
%
%  Required inputs
%  ---------------
%  mesh_mat :
%      Path to the tetrahedral head mesh MAT file. Used to load node and
%      element arrays and determine the set of brain nodes.
%
%  nirsmodel_file :
%      Path to the forward-model MAT file. Expected to contain nirs_mesh and
%      wavelength-specific forward-model components required by
%      fnirspower.fwdcomp.buildPathAdjJacobian.
%
%  roi_def :
%      Struct defining one or more ROI seeds. Required field:
%
%        roi_def.centers_mm :
%            ROI center coordinates in millimeters, either:
%              [1 x 3]   for a single ROI
%              [nROI x 3] for multiple ROIs
%
%      Optional fields:
%
%        roi_def.radius_mm :
%            ROI radius/radii used for sphere-based ROI seeds.
%
%        roi_def.sigma_mm :
%            ROI sigma/sigmas used for Gaussian ROI seeds and subject-level
%            ROI generation.
%
%        roi_def.amp_std :
%            ROI fractional amplitude standard deviation used when 
%            converting ROI maps to Hb/Δμa maps.
%
%        roi_def.names :
%            Cell array of ROI names.
%
%  dmua_file :
%      Path to a saved absorption-estimation result file containing
%      a spectral absorption magnitude estimate. Pass '' to omit and provide
%      opts.dmu_a directly instead.
%
%  Name-value options
%  ------------------
%  Core simulation options
%
%  'nSubjects' :
%      Number of simulated subjects used when generating ROI probability
%      maps. Default: 80
%
%  'hbRatio' :
%      HbO:HbR ratio used in ROI-to-Hb conversion. Default: 3
%
%  'dmuaScale' :
%      If true, divide the resolved dmu_a magnitude by opts.dmuaScaleDiv
%      before prediction. Default: true
%
%  'dmuaScaleDiv' :
%      Divisor applied to dmu_a when opts.dmuaScale is true. Default: 3
%
%  'epsilon_mm_uM' :
%      Extinction coefficient matrix used to build the spectral Jacobian.
%      Default:
%        2.303e-6 * [58.6,154.8; 105.8,69.1]
%
%  'dmu_a' :
%      Direct spectral absorption magnitude vector:
%        [dmu760; dmu850]
%      Used if dmua_file is not provided, or if opts.dmu_a_override is true.
%      Default: []
%
%  'dmu_a_override' :
%      If true, opts.dmu_a takes precedence over dmua_file even when both
%      are supplied. Default: false
%
%  ROI/mesh visualization options
%
%  'plotMesh' :
%      If true, plot the brain mesh. Default: false
%
%  'plotROISeed' :
%      If true, plot the template ROI seed(s). Default: true
%
%  'plotAvgProb' :
%      If true, plot the average probabilistic ROI map for each ROI.
%      Default: true
%
%  'plotSubjects' :
%      Number of example subject-specific ROI maps to plot. Default: 0
%
%  'roiSeedMode' :
%      ROI seed type for visualization. Supported values:
%        'sphere'
%        'gaussian'
%      Default: 'sphere'
%
%  'cutoff_axis' :
%      Axis used for mesh cutoff visualization. Default: 'y'
%
%  'cutoff_value' :
%      Cutoff value used for mesh visualization. Default: 100
%
%  'edge_alpha' :
%      Edge transparency for mesh plotting. Default: 0.3
%
%  'face_alpha' :
%      Face transparency for mesh plotting. Default: 1
%
%  'view_angle' :
%      Two-element viewing angle for mesh plots. Default: [180 90]
%
%  Channel-plotting and layout options
%
%  'plotChannels' :
%      If true, plot mean predicted channel-level HbO/HbR maps for each ROI.
%      Default: false
%
%  'montageName' :
%      Montage name used to help resolve a layout file when plotChannels is
%      true and opts.layout_file is not supplied. Default: ''
%
%  'layout_file' :
%      Explicit path to a layout MAT file used for channel-level plotting.
%      Default: ''
%
%  'layout_search_dir' :
%      Directory searched for a matching layout file when layout_file is not
%      provided. Default: ''
%
%  'layout_var' :
%      Optional variable name to load from layout_file. If empty, the first
%      variable in the MAT file is used. Default: ''
%
%  'channel_clim' :
%      Optional color-axis limit for channel-level plots. If empty, the
%      limit is determined from the data. Default: []
%
%  Distance computation
%
%  'computeDistances' :
%      If true, compute distances between channel midpoints and ROI
%      centroids. Default: true
%
%  Returns
%  -------
%  M : struct with fields
%
%    M.roi :
%        Struct containing ROI definitions and simulated ROI maps:
%
%          .names :
%              ROI names
%
%          .centers_mm :
%              ROI center coordinates
%
%          .radius_mm :
%              ROI radii
%
%          .sigma_mm :
%              ROI sigmas
%
%          .amp_std :
%              ROI amplitude variability
%
%          .avg_prob :
%              Average ROI probability map for each ROI [nROI x nBrain]
%
%          .subj_prob :
%              Subject-level ROI probability maps
%              [nSubjects x nROI x nBrain]
%
%          .centroids_mm :
%              Probability-weighted ROI centroids [nROI x 3]
%
%          .info :
%              Cell array of ROI-generation metadata returned by
%              fnirspower.measpred.compute_subject_ROI
%
%    M.channels :
%        Struct containing channel-space geometry:
%
%          .pos_mm :
%              Channel midpoint coordinates [nCh x 3]
%
%          .dist2roi_mm :
%              Distance from each channel midpoint to each ROI centroid
%              [nCh x nROI]
%
%    M.vox :
%        Struct containing voxel/brain-node-level predicted quantities:
%
%          .HbO
%          .HbR
%          .muA760
%          .muA850
%
%    M.chan :
%        Struct containing channel-level predicted hemoglobin responses:
%
%          .HbO
%          .HbR
%
%    M.meta :
%        Struct containing metadata and resolved inputs:
%
%          .mesh_mat
%          .nirsmodel_file
%          .dmua_file
%          .layout_file
%          .dmu_a_input
%          .dmu_a_used
%          .opts
%          .SJ
%
%  Example
%  -------
%  roi_def = struct();
%  roi_def.centers_mm = [-70 -10 9];
%  roi_def.radius_mm = 2;
%  roi_def.sigma_mm = 4;
%  roi_def.amp_std = 0.3;
%  roi_def.names = {'AuditoryDepth'};
%
%  M = fnirspower.pipeline.run_measurement_prediction( ...
%      mesh_mat, nirsmodel_file, roi_def, '', ...
%      'dmu_a', [1e-3; 1e-3], ...
%      'dmuaScale', false, ...
%      'nSubjects', 80, ...
%      'hbRatio', 3, ...
%      'plotMesh', true, ...
%      'plotROISeed', true, ...
%      'plotAvgProb', true, ...
%      'plotChannels', true, ...
%      'layout_file', layout_file);
%
%  Example using a saved absorption-estimation result
%  --------------------------------------------------
%  M = fnirspower.pipeline.run_measurement_prediction( ...
%      mesh_mat, nirsmodel_file, roi_def, dmua_file, ...
%      'nSubjects', 80, ...
%      'hbRatio', 3, ...
%      'plotChannels', true, ...
%      'layout_file', layout_file);
%
%  Notes
%  -----
%  - If dmua_file is empty, opts.dmu_a must be provided unless
%    opts.dmu_a_override is false and a valid file is supplied.
%  - Layout resolution for channel plotting proceeds in this order:
%      1) opts.layout_file
%      2) opts.layout_search_dir
%      3) alongside nirsmodel_file
%      4) fnirspower.paths().layouts (fallback)
%  - The function does not assume a layout-specific variable name inside the
%    layout MAT file; if opts.layout_var is empty, the first variable is used.
%  - Channel-level predictions are noiseless forward projections based on the
%    supplied or estimated spectral absorption magnitude and simulated ROI
%    distributions.
%
%  See also
%  --------
%  fnirspower.measpred.compute_subject_ROI
%  fnirspower.measpred.ROImua_to_HbOchan
%  fnirspower.fwdcomp.buildPathAdjJacobian
%  fnirspower.measpred.plot_brain_mesh_iso2mesh
%  fnirspower.helpers.plot_channel_level_hb

arguments
    mesh_mat (1,:) char
    nirsmodel_file (1,:) char
    roi_def struct
    dmua_file (1,:) char = ''
    
    opts.nSubjects (1,1) double = 80
    opts.hbRatio (1,1) double = 3
    opts.dmuaScale (1,1) logical = true
    opts.dmuaScaleDiv (1,1) double = 3
    opts.epsilon_mm_uM double = 2.303e-6 * [58.6,154.8; 105.8,69.1]
    
    opts.dmu_a double = []
    opts.dmu_a_override (1,1) logical = false
    
    opts.plotMesh (1,1) logical = false
    opts.plotROISeed (1,1) logical = true
    opts.plotAvgProb (1,1) logical = true
    opts.plotSubjects (1,1) double = 0
    opts.roiSeedMode (1,:) char = 'sphere'
    opts.cutoff_axis (1,:) char = 'y'
    opts.cutoff_value (1,1) double = 100
    opts.edge_alpha (1,1) double = 0.3
    opts.face_alpha (1,1) double = 1
    opts.view_angle (1,2) double = [180 90]
    
    opts.plotChannels (1,1) logical = false
    opts.montageName (1,:) char = ''
    opts.layout_file (1,:) char = ''
    opts.layout_search_dir (1,:) char = ''
    opts.layout_var (1,:) char = ''
    opts.channel_clim double = []
    
    opts.computeDistances (1,1) logical = true
end

layout = [];
if opts.plotChannels
    if isempty(opts.layout_file)
        if isempty(opts.montageName)
            error(['plotChannels=true requires either opts.layout_file or ', ...
                'opts.montageName (with an accessible search location).']);
        end
        opts.layout_file = local_resolve_layout_file(opts.montageName, opts.layout_search_dir, nirsmodel_file);
    end
    
    L = load(opts.layout_file);
    if ~isempty(opts.layout_var)
        if ~isfield(L, opts.layout_var)
            error('layout_var="%s" not found in %s', opts.layout_var, opts.layout_file);
        end
        layout = L.(opts.layout_var);
    else
        f = fieldnames(L);
        if isempty(f)
            error('No variables found in layout file: %s', opts.layout_file);
        end
        layout = L.(f{1});
        opts.layout_var = f{1};
    end
end

% ---- load mesh + brain indices
mesh_model = load(mesh_mat);
if isfield(mesh_model,'genmesh')
    p = mesh_model.genmesh.node;
    e = mesh_model.genmesh.ele;
else
    p = mesh_model.p;
    e = mesh_model.e;
end

brain_tet = e(e(:,5)==2,:);
brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);
p_brain = p(brain_nodes_idx,:);
nBrain  = numel(brain_nodes_idx);

% ---- load model + channel positions
nirs_model = load(nirsmodel_file, 'nirs_mesh', 'J_760', 'J_850', 'DATA760', 'DATA850');
if ~isfield(nirs_model,'nirs_mesh')
    error('nirsmodel_file must contain nirs_mesh.');
end

link_source   = nirs_model.nirs_mesh.link(:,1);
link_detector = nirs_model.nirs_mesh.link(:,2);
src_pos = nirs_model.nirs_mesh.source.coord;
det_pos = nirs_model.nirs_mesh.meas.coord;

nCh = numel(link_source);
channel_pos = zeros(nCh,3);
for ch = 1:nCh
    channel_pos(ch,:) = (src_pos(link_source(ch),:) + det_pos(link_detector(ch),:))/2;
end

% ---- resolve spectral amplitude dmu_a
dmua_path = '';
dmu_a = [];

if opts.dmu_a_override
    dmu_a = opts.dmu_a(:);
end

if isempty(dmu_a) && ~isempty(dmua_file)
    dmua_path = char(dmua_file);
    if exist(dmua_path,'file') ~= 2
        error('Specified dmua_file does not exist: %s', dmua_path);
    end
    
    D = load(dmua_path);
    if isfield(D,'dmu_a_final')
        dmu_a = D.dmu_a_final;
    end
    if isempty(dmu_a) && isfield(D,'dmu_a_opt')
        dmu_a = D.dmu_a_opt;
    end
    if isempty(dmu_a) && isfield(D,'out') && isfield(D.out,'dmu_a_opt')
        dmu_a = D.out.dmu_a_opt;
    end
    if isempty(dmu_a)
        error('Could not find dmu_a in dmua_file: %s', dmua_path);
    end
end

if isempty(dmu_a) && ~isempty(opts.dmu_a)
    dmu_a = opts.dmu_a(:);
end

if isempty(dmu_a)
    error(['No dmu_a source available. Provide either dmua_file or ', ...
        '''dmu_a'', [dmu760; dmu850].']);
end

if numel(dmu_a) ~= 2
    error('dmu_a must be 2x1 [760; 850].');
end

if opts.dmuaScale
    dmu_a_used = dmu_a / opts.dmuaScaleDiv;
else
    dmu_a_used = dmu_a;
end

% ---- ROI definition
if ~isfield(roi_def,'centers_mm') || isempty(roi_def.centers_mm)
    error('roi_def.centers_mm is required (1x3 or nROI x 3).');
end

centers = roi_def.centers_mm;
if isvector(centers) && numel(centers)==3
    centers = reshape(centers,1,3);
end
if size(centers,2)~=3
    error('roi_def.centers_mm must be [nROI x 3].');
end
nROI = size(centers,1);

radius = get_roi_vec(roi_def, 'radius_mm', nROI, 7.5);
radius_std = get_roi_vec(roi_def, 'radius_std_mm', nROI, 0);
sigma  = get_roi_vec(roi_def, 'sigma_mm',  nROI, 5);
amp_std    = get_roi_vec(roi_def, 'amp_std',       nROI, 0.3);

names = {};
if isfield(roi_def,'names')
    names = roi_def.names;
end
if isempty(names)
    names = arrayfun(@(k) sprintf('ROI_%d',k), 1:nROI, 'UniformOutput', false);
end

% ---- generate subject ROI probability maps on brain nodes (PER ROI)
subj_prob = NaN(opts.nSubjects, nROI, nBrain);
avg_prob  = NaN(nROI, nBrain);
roi_info  = cell(nROI,1);

for r = 1:nROI
    def_r = struct();

    def_r.center_mm = centers(r, :);
    def_r.radius_mm = radius(r);
    def_r.radius_std_mm = radius_std(r);
    def_r.sigma_mm = sigma(r);

    roi_opts_r = struct();
    roi_opts_r.n_subj = opts.nSubjects;

    [subj_r, avg_r, info_r] = ...
        fnirspower.measpred.compute_subject_ROI( ...
            p_brain, ...
            def_r, ...
            roi_opts_r);

    subj_prob(:, r, :) = subj_r;
    avg_prob(r, :) = avg_r(:)';
    roi_info{r} = info_r;
end

% ---- OPTIONAL iso2mesh plots
plot_cfg = struct( ...
    'cutoff_axis',  opts.cutoff_axis, ...
    'cutoff_value', opts.cutoff_value, ...
    'edge_alpha',   opts.edge_alpha, ...
    'face_alpha',   opts.face_alpha, ...
    'view_angle',   opts.view_angle);

if opts.plotMesh
    fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
        p, e, brain_nodes_idx, [], plot_cfg, ...
        'title_str', 'Brain mesh');
end

if opts.plotROISeed
    for r = 1:nROI
        seed = zeros(nBrain,1);
        
        if strcmpi(opts.roiSeedMode,'sphere')
            d = sqrt(sum((p_brain - centers(r,:)).^2, 2));
            seed(d <= radius(r)) = 1;
        elseif strcmpi(opts.roiSeedMode,'gaussian')
            d2 = sum((p_brain - centers(r,:)).^2, 2);
            seed = exp(-d2/(2*sigma(r)^2));
            seed = seed / max(eps, max(seed));
        else
            error('opts.roiSeedMode must be "sphere" or "gaussian".');
        end
        
        fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
            p, e, brain_nodes_idx, seed, plot_cfg, ...
            'title_str', sprintf('ROI seed (%s): %s', opts.roiSeedMode, names{r}));
    end
end

if opts.plotAvgProb
    for r = 1:nROI
        fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
            p, e, brain_nodes_idx, avg_prob(r,:).', plot_cfg, ...
            'title_str', sprintf('Avg ROI probability: %s', names{r}));
    end
end

if opts.plotSubjects > 0
    nShow = min(opts.plotSubjects, opts.nSubjects);
    for s = 1:nShow
        for r = 1:nROI
            fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
                p, e, brain_nodes_idx, squeeze(subj_prob(s,r,:)), plot_cfg, ...
                'title_str', sprintf('Subject %d ROI prob: %s', s, names{r}));
        end
    end
end

% ---- spectral Jacobian for forward projection
SJ = fnirspower.fwdcomp.buildPathAdjJacobian(nirs_model, opts.epsilon_mm_uM, brain_nodes_idx);

% ---- per subject + ROI forward prediction
vox_HbO   = NaN(opts.nSubjects, nROI, nBrain);
vox_HbR   = NaN(opts.nSubjects, nROI, nBrain);
vox_mu760 = NaN(opts.nSubjects, nROI, nBrain);
vox_mu850 = NaN(opts.nSubjects, nROI, nBrain);

chan_HbO  = NaN(opts.nSubjects, nROI, nCh);
chan_HbR  = NaN(opts.nSubjects, nROI, nCh);

amplitude_scales = NaN(opts.nSubjects, nROI);

for s = 1:opts.nSubjects
    for r = 1:nROI
        roi_prob_brain = squeeze(subj_prob(s,r,:));
        
        [HbO_full, HbR_full, muA760_full, muA850_full, ...
            amplitude_scales(s, r)] = ...
            fnirspower.measpred.ROImua_to_HbOchan( ...
            p, e, brain_nodes_idx, dmu_a_used, opts.hbRatio, ...
            'roi_prob', roi_prob_brain, ...
            'roi_amp_std', amp_std(r));
        
        vox_HbO(s,r,:)   = HbO_full(brain_nodes_idx);
        vox_HbR(s,r,:)   = HbR_full(brain_nodes_idx);
        vox_mu760(s,r,:) = muA760_full(brain_nodes_idx);
        vox_mu850(s,r,:) = muA850_full(brain_nodes_idx);
        
        mu_vec = [muA760_full(brain_nodes_idx); muA850_full(brain_nodes_idx)];
        dHb_ch = SJ.J_chrom_brain * mu_vec;
        chan_HbO(s,r,:) = dHb_ch(1:nCh);
        chan_HbR(s,r,:) = dHb_ch(nCh+1:end);
    end
end

% ---- Channel-level plots
if opts.plotChannels
    if isempty(layout)
        error('Internal error: layout could not be resolved.');
    end
    
    for r = 1:nROI
        HbO_mean = squeeze(mean(chan_HbO(:,r,:), 1, 'omitnan'));
        HbR_mean = squeeze(mean(chan_HbR(:,r,:), 1, 'omitnan'));
        data2plot = [HbO_mean(:)'; HbR_mean(:)'];
        
        if isempty(opts.channel_clim)
            clim = max(abs(data2plot(:)));
        else
            clim = opts.channel_clim;
        end
        
        fnirspower.helpers.plot_channel_level_hb( ...
            layout, data2plot, clim, datestr(datetime('now'),'yyyy-mm-dd'), ...
            sprintf('Predicted-Hb-%s', names{r}), 'MeasurementPrediction', 1);
    end
end

% ---- distances
ROI_centroids = NaN(nROI,3);
dist2roi      = NaN(nCh,nROI);

if opts.computeDistances
    for r = 1:nROI
        pr = avg_prob(r,:).';
        pr = pr / max(eps, sum(pr));
        ROI_centroids(r,:) = (pr.' * p_brain);
        diffs = channel_pos - ROI_centroids(r,:);
        dist2roi(:,r) = sqrt(sum(diffs.^2,2));
    end
end

% ---- pack outputs
M = struct();

M.roi = struct();
M.roi.names        = names;
M.roi.centers_mm   = centers;
M.roi.radius_mm    = radius;
M.roi.radius_std_mm = radius_std;
M.roi.sigma_mm     = sigma;
M.roi.amp_std          = amp_std;
M.roi.amplitude_scales = amplitude_scales;
M.roi.avg_prob     = avg_prob;
M.roi.subj_prob    = subj_prob;
M.roi.centroids_mm = ROI_centroids;
M.roi.info         = roi_info;

M.channels = struct();
M.channels.pos_mm      = channel_pos;
M.channels.dist2roi_mm = dist2roi;

M.vox  = struct('HbO',vox_HbO,'HbR',vox_HbR,'muA760',vox_mu760,'muA850',vox_mu850);
M.chan = struct('HbO',chan_HbO,'HbR',chan_HbR);

M.meta = struct();
M.meta.mesh_mat       = char(mesh_mat);
M.meta.nirsmodel_file = char(nirsmodel_file);
M.meta.dmua_file     = dmua_path;
M.meta.layout_file    = char(opts.layout_file);
M.meta.dmu_a_input    = dmu_a;
M.meta.dmu_a_used     = dmu_a_used;
M.meta.opts           = opts;

end


function layout_file = local_resolve_layout_file(montageName, layout_search_dir, nirsmodel_file)

cand_names = { ...
    sprintf('layout_%s.mat', montageName), ...
    sprintf('layout_*%s*.mat', montageName)};

search_dirs = {};

if ~isempty(layout_search_dir)
    search_dirs{end+1} = layout_search_dir;
end

nirsmodel_dir = fileparts(nirsmodel_file);
if ~isempty(nirsmodel_dir)
    search_dirs{end+1} = nirsmodel_dir;
end

% Final fallback to package paths
try
    P = fnirspower.paths();
    if isfield(P,'layouts') && ~isempty(P.layouts)
        search_dirs{end+1} = P.layouts;
    end
catch
end

for di = 1:numel(search_dirs)
    this_dir = search_dirs{di};
    if isempty(this_dir) || exist(this_dir,'dir')~=7
        continue;
    end
    
    exact_file = fullfile(this_dir, cand_names{1});
    if exist(exact_file,'file') == 2
        layout_file = exact_file;
        return
    end
    
    d = dir(fullfile(this_dir, cand_names{2}));
    if ~isempty(d)
        [~,ix] = sort([d.datenum], 'descend');
        layout_file = fullfile(d(ix(1)).folder, d(ix(1)).name);
        return
    end
end

error(['Could not resolve layout file for montageName="%s". Provide opts.layout_file ', ...
    'or opts.layout_search_dir.'], montageName);

end


function v = get_roi_vec(roi_def, field, nROI, defaultVal)
if isfield(roi_def, field) && ~isempty(roi_def.(field))
    v = roi_def.(field);
else
    v = defaultVal;
end

if isscalar(v)
    v = repmat(v, nROI, 1);
else
    v = v(:);
    if numel(v) ~= nROI
        error('roi_def.%s must be scalar or length nROI (%d).', field, nROI);
    end
end
end