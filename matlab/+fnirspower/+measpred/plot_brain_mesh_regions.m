function hFig = plot_brain_mesh_regions(p, e, brain_nodes_idx, roi_maps, roi_centers_mm, roi_names, varargin)
%FNIRSPOWER.MEASPRED.PLOT_BRAIN_MESH_REGIONS
% Plot brain mesh (iso2mesh) and optionally overlay ROI maps + ROI centers.
%
% Required
%   p              : nNodes x 3 node coordinates
%   e              : nElem x 5 tetrahedra (col 5 = region id)
%   brain_nodes_idx: nBrain x 1 indices of brain nodes in p
%   roi_maps       : nROI x nBrain (or nBrain x 1) ROI probability/weight maps on brain nodes
%   roi_centers_mm : nROI x 3 (can be [])
%   roi_names      : 1 x nROI cellstr/string (can be [])
%
% Name-value options
%   'region_id'      (default 2)
%   'cutoff_axis'    (default 'y')
%   'cutoff_value'   (default 100)
%   'edge_alpha'     (default 0.3)
%   'face_alpha'     (default 1)
%   'view_angle'     (default [180 90])
%   'figure_pos'     (default [600 200 1200 800])
%   'brain_face_rgb' (default [] -> brewermap(9,'Spectral')(3,:) if available)
%   'plot_roi'       (default true)
%   'roi_threshold'  (default 0.25)  % fraction of max per ROI used to show ROI nodes
%   'roi_point_size' (default 12)
%   'plot_centers'   (default true)
%   'center_size'    (default 80)
%
% Returns
%   hFig : figure handle

% -----------------------
% input parsing
% -----------------------
ip = inputParser;
ip.FunctionName = mfilename;

addRequired(ip,'p',@(x)isnumeric(x) && size(x,2)==3);
addRequired(ip,'e',@(x)isnumeric(x) && size(x,2)>=5);
addRequired(ip,'brain_nodes_idx',@(x)isnumeric(x) && isvector(x));
addRequired(ip,'roi_maps',@(x)isnumeric(x) && ~isempty(x));
addRequired(ip,'roi_centers_mm',@(x)isnumeric(x) || isempty(x));
addRequired(ip,'roi_names',@(x)iscell(x) || isstring(x) || isempty(x));

addParameter(ip,'region_id',2,@(x)isnumeric(x) && isscalar(x));
addParameter(ip,'cutoff_axis','y',@(x)ischar(x) || (isstring(x)&&isscalar(x)));
addParameter(ip,'cutoff_value',100,@(x)isnumeric(x) && isscalar(x));
addParameter(ip,'edge_alpha',0.3,@(x)isnumeric(x) && isscalar(x));
addParameter(ip,'face_alpha',1,@(x)isnumeric(x) && isscalar(x));
addParameter(ip,'view_angle',[180 90],@(x)isnumeric(x) && numel(x)==2);
addParameter(ip,'figure_pos',[600 200 1200 800],@(x)isnumeric(x) && numel(x)==4);
addParameter(ip,'brain_face_rgb',[],@(x)isnumeric(x) && (isempty(x) || numel(x)==3));
addParameter(ip,'plot_roi',true,@(x)islogical(x) && isscalar(x));
addParameter(ip,'roi_threshold',0.25,@(x)isnumeric(x) && isscalar(x) && x>=0 && x<=1);
addParameter(ip,'roi_point_size',12,@(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(ip,'plot_centers',true,@(x)islogical(x) && isscalar(x));
addParameter(ip,'center_size',80,@(x)isnumeric(x) && isscalar(x) && x>0);

parse(ip,p,e,brain_nodes_idx,roi_maps,roi_centers_mm,roi_names,varargin{:});
O = ip.Results;

% -----------------------
% normalize ROI inputs
% -----------------------
brain_nodes_idx = brain_nodes_idx(:);

if isvector(roi_maps)
  roi_maps = roi_maps(:).'; % 1 x nBrain
end
if size(roi_maps,2) ~= numel(brain_nodes_idx)
  error('roi_maps must be nROI x nBrain (brain-only). Got %d cols, expected %d.', size(roi_maps,2), numel(brain_nodes_idx));
end
nROI = size(roi_maps,1);

if isempty(roi_names)
  roi_names = arrayfun(@(k)sprintf('ROI_%d',k), 1:nROI, 'UniformOutput', false);
elseif isstring(roi_names)
  roi_names = cellstr(roi_names);
end
if numel(roi_names) ~= nROI
  error('roi_names must have length nROI (%d).', nROI);
end

if ~isempty(roi_centers_mm)
  if size(roi_centers_mm,2)~=3 || size(roi_centers_mm,1)~=nROI
    error('roi_centers_mm must be nROI x 3.');
  end
end

% -----------------------
% colors
% -----------------------
brain_rgb = O.brain_face_rgb;
if isempty(brain_rgb)
  try
    map = brewermap(9,'Spectral');
    brain_rgb = map(3,:);
  catch
    brain_rgb = [0.2 0.6 0.6]; % fallback
  end
end

% ROI colors (distinct)
roi_colors = lines(max(nROI,1));

% cutoff string like: 'y<100'
cutoff_ax  = char(O.cutoff_axis);
cutoff_txt = sprintf('%s<%g', cutoff_ax, O.cutoff_value);

% -----------------------
% plot
% -----------------------
hFig = figure('Name','Brain mesh + ROI overlay');
set(hFig,'Position',O.figure_pos);
hold on;

% brain tets
tet = e(e(:,5)==O.region_id, :);
iso2mesh_plotmesh(p, tet, cutoff_txt, ...
  'FaceColor', brain_rgb, 'EdgeAlpha', O.edge_alpha, 'FaceAlpha', O.face_alpha);

xlabel('x'); ylabel('y'); zlabel('z');
view(O.view_angle);
axis tight;
axis vis3d;

% overlay ROI nodes as points
if O.plot_roi
  brain_xyz = p(brain_nodes_idx,:);
  for r = 1:nROI
    w = roi_maps(r,:).';
    wmax = max(w);
    if ~isfinite(wmax) || wmax<=0
      continue;
    end
    thr = O.roi_threshold * wmax;
    keep = w >= thr;

    scatter3(brain_xyz(keep,1), brain_xyz(keep,2), brain_xyz(keep,3), ...
      O.roi_point_size, roi_colors(r,:), 'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeAlpha', 0.0);
  end
end

% overlay ROI centers
if O.plot_centers && ~isempty(roi_centers_mm)
  for r = 1:nROI
    scatter3(roi_centers_mm(r,1), roi_centers_mm(r,2), roi_centers_mm(r,3), ...
      O.center_size, roi_colors(r,:), 'p', 'filled', 'MarkerEdgeColor','k');
  end
end

% legend
leg = [{'Brain'}, roi_names(:)'];
legend(leg{:}, 'Location','bestoutside');

hold off;
end