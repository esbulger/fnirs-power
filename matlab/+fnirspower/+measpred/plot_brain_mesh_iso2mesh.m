function plot_brain_mesh_iso2mesh(p, e, brain_nodes_idx, overlay_brain, plot_cfg, varargin)
%PLOT_BRAIN_MESH_ISO2MESH Plot a brain mesh with optional nodewise overlay using iso2mesh.
%
%  fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
%      p, e, brain_nodes_idx, overlay_brain, plot_cfg)
%
%  fnirspower.measpred.plot_brain_mesh_iso2mesh( ...
%      p, e, brain_nodes_idx, overlay_brain, plot_cfg, ...
%      'title_str', title_str)
%
%  Description
%  -----------
%  This helper plots the brain portion of a tetrahedral head mesh using
%  iso2mesh_plotmesh. If a brain-node scalar field is provided, the values are
%  embedded into the full node array and displayed on the plotted brain mesh.
%
%  Inputs
%  ------
%  p :
%      Full mesh node array [nNodes x 3].
%
%  e :
%      Full mesh element array [nElem x >=5]. Column 5 is assumed to contain
%      tissue/region labels, with label 2 corresponding to brain.
%
%  brain_nodes_idx :
%      Indices of brain nodes within p.
%
%  overlay_brain :
%      Either [] for mesh-only plotting, or a brain-only scalar vector
%      [nBrain x 1] to be overlaid on the brain mesh.
%
%  plot_cfg :
%      Struct containing plotting options. Expected fields:
%        .cutoff_axis   - axis for iso2mesh cutoff string, e.g. 'y'
%        .cutoff_value  - cutoff value, e.g. 100
%        .edge_alpha    - edge transparency
%        .face_alpha    - face transparency
%        .view_angle    - 1x2 view angle, e.g. [180 90]
%
%  Name-value inputs
%  -----------------
%  'title_str' :
%      Optional figure title. Default: ''
%
%  Notes
%  -----
%  - Brain elements are taken as e(e(:,5)==2,:).
%  - If overlay_brain is provided, the color axis is centered at zero using
%    symmetric limits [-max(abs(v)) max(abs(v))].
%  - The function uses the redblue colormap if it is on the MATLAB path.
%
%  See also
%  --------
%  iso2mesh_plotmesh
%  redblue

ip = inputParser;
addParameter(ip, 'title_str', '', @(x)ischar(x) || isstring(x));
parse(ip, varargin{:});
title_str = char(ip.Results.title_str);

% Cutoff string for iso2mesh
cutoff_txt = sprintf('%s<%g', plot_cfg.cutoff_axis, plot_cfg.cutoff_value);

% Brain elements only
tet_gm = e(e(:,5)==2, :);
tet_wm = e(e(:,5)==1, :);

% Prepare plotting array
if isempty(overlay_brain)
    plot_points = p;
else
    overlay_brain = overlay_brain(:);
    if numel(overlay_brain) ~= numel(brain_nodes_idx)
        error('overlay_brain must be length nBrain (%d).', numel(brain_nodes_idx));
    end
    v = zeros(size(p,1),1);
    v(brain_nodes_idx) = overlay_brain;
    plot_points = [p, v];
end

% Create figure
figure;
set(gcf, 'Position', [600 200 1200 800]);
hold on

% Base face color for mesh-only plotting
face_col = [0.9 0.9 0.9];

iso2mesh_plotmesh(plot_points, tet_gm, cutoff_txt, ...
    'EdgeAlpha', plot_cfg.edge_alpha, ...
    'FaceAlpha', plot_cfg.face_alpha);
iso2mesh_plotmesh(p, tet_wm, cutoff_txt, ...
    'FaceColor', face_col, ...
    'EdgeAlpha', plot_cfg.edge_alpha, ...
    'FaceAlpha', plot_cfg.face_alpha);
xlabel('x');
ylabel('y');
zlabel('z');
view(plot_cfg.view_angle(1), plot_cfg.view_angle(2));
axis equal
axis off

% Overlay settings
if ~isempty(overlay_brain)
    if exist('redblue', 'file') == 2
        colormap(redblue);
    end
    colorbar;

    clim = max(abs(overlay_brain(:)));
    if isfinite(clim) && clim > 0
        caxis([-clim, clim]);
    end
end

if ~isempty(title_str)
    title(title_str, 'Interpreter', 'none');
end
end