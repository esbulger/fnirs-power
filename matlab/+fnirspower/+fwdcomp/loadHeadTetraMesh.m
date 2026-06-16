function head = loadHeadTetraMesh(meshMatPath)
%LOADHEADTETRAMESH Load a 5-layer head tetrahedral mesh produced upstream.
% Returns struct with fields:
% node (Nx3), ele (Mx5+) where column 5 is region id (1..5)
% Notes
% - This wraps whatever upstream saved as genmesh/node/ele into a standard form.

S = load(meshMatPath);

if isfield(S,'genmesh')
head.node = S.genmesh.node;
head.ele = S.genmesh.ele;
elseif all(isfield(S, {'node','elem'}))
head.node = S.node; % iso2mesh style
head.ele = S.elem;
else
error('Unrecognized mesh file format: %s', meshMatPath);
end

if size(head.ele,2) < 5
error('Expected region labels in column 5 of head.ele.');
end
end