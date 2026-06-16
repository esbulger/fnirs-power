function [srcProj, detProj, link] = placeAndProjectOptodes(head, probeInfo, opts)
%PLACEANDPROJECTOPTODES Adjust probe coords and project to scalp nodes.
%
%  [srcProj, detProj, link] = fnirspower.fwdcomp.placeAndProjectOptodes(head, probeInfo, opts)
%
%  Inputs
%  ------
%  head      : struct with fields
%                .node  [N x 3]
%                .ele   [M x >=5], where column 5 is tissue/region label
%
%  probeInfo : struct with fields
%                .probes.coords_s3   [S x 3]
%                .probes.coords_d3   [D x 3]
%                .probes.index_c     [C x 2]
%
%  opts      : struct with optional fields
%                .scale        (default 10)
%                .offset       (default [0 18 0])
%                .zero_center  (default false)
%                .debugPlots   (default false)
%                .cutoff_axis  (default 'y')
%                .cutoff_value (default 20)
%                .edge_alpha   (default 0.2)
%                .face_alpha   (default 1)
%
%  Outputs
%  -------
%  srcProj   : projected source coordinates on scalp nodes
%  detProj   : projected detector coordinates on scalp nodes
%  link      : channel link matrix copied from probeInfo.probes.index_c

arguments
    head struct
    probeInfo struct
    opts struct = struct()
end

% Defaults
if ~isfield(opts, 'scale'),        opts.scale = 10; end
if ~isfield(opts, 'offset'),       opts.offset = [0 18 0]; end
if ~isfield(opts, 'zero_center'),  opts.zero_center = false; end
if ~isfield(opts, 'debugPlots'),   opts.debugPlots = false; end
if ~isfield(opts, 'cutoff_axis'),  opts.cutoff_axis = 'y'; end
if ~isfield(opts, 'cutoff_value'), opts.cutoff_value = 150; end
if ~isfield(opts, 'edge_alpha'),   opts.edge_alpha = 0.2; end
if ~isfield(opts, 'face_alpha'),   opts.face_alpha = 0.9; end

p = head.node;
e = head.ele;

src = probeInfo.probes.coords_s3;
det = probeInfo.probes.coords_d3;
opt = [src; det];

if opts.zero_center
    mu = mean(opt, 1);
else
    mu = [0 0 0];
end

srcRot = opts.scale * src - mu + opts.offset;
detRot = opts.scale * det - mu + opts.offset;

% Use outer scalp nodes only for projection.
scalpTet = e(e(:,5)==5, :);
scalpNodeIdx = unique(scalpTet(:,1:4));
scalpNodes = p(scalpNodeIdx, :);

% Project to closest scalp node
srcIdxLocal = dsearchn(scalpNodes, srcRot);
detIdxLocal = dsearchn(scalpNodes, detRot);

srcProj = scalpNodes(srcIdxLocal, :);
detProj = scalpNodes(detIdxLocal, :);

if opts.debugPlots
    cutoff_txt = [opts.cutoff_axis '<' num2str(opts.cutoff_value)];

    map = [
        0.70, 0.50, 0.90;   % region 4
        0.40, 0.20, 0.60    % region 5
    ];

    figure('Name','Optode placement debug');
    set(gcf,'Position',[600 200 1200 800]);
    hold on

    if any(e(:,5)==5)
        iso2mesh_plotmesh(p, e(e(:,5)==5,:), cutoff_txt, ...
            'FaceColor', map(2,:), ...
            'EdgeAlpha', opts.edge_alpha, ...
            'FaceAlpha', opts.face_alpha);
    end

    % Plot raw transformed optodes
    scatter3(srcRot(:,1), srcRot(:,2), srcRot(:,3), 60, 'r', 'filled');
    scatter3(detRot(:,1), detRot(:,2), detRot(:,3), 60, 'b', 'filled');

    % Plot projected scalp-node positions
    scatter3(srcProj(:,1), srcProj(:,2), srcProj(:,3), 70, 'k', 'filled');
    scatter3(detProj(:,1), detProj(:,2), detProj(:,3), 70, 'm', 'filled');

    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    view(-90, 0)
    legend({'Region 4','Region 5','src raw','det raw','src proj','det proj'});
end

link = probeInfo.probes.index_c;
end