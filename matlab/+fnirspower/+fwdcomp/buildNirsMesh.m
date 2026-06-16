function nirs = buildNirsMesh(baseMeshPath, srcCoord, detCoord, link, outPath, opts)
%BUILDNIRSMESH Load a NIRFASTer mesh, place optodes, and save.
% baseMeshPath : path to a NIRFASTer .mat mesh (for load_mesh)
% srcCoord : Sx3 double, source coordinates already projected to scalp
% detCoord : Dx3 double, detector coordinates already projected to scalp
% link : Cx2 channel matrix [srcIdx detIdx]
% outPath : path (w/o extension) to save the updated mesh via save_mesh
% opts : struct with fields (optional):
% ri : refractive index scalar or per-node vector (default 1.4)
% fwhm : Sx1 source FWHM (mm), default zeros
%
arguments
    baseMeshPath (1,:) char
    srcCoord double
    detCoord double
    link double
    outPath (1,:) char
    opts struct = struct()
end

if ~isfield(opts, 'ri') || isempty(opts.ri)
    opts.ri = 1.4;
end
if ~isfield(opts, 'fwhm') || isempty(opts.fwhm)
    opts.fwhm = [];
end

nirs = load_mesh(baseMeshPath);
nirs.source.coord = srcCoord;
nirs.meas.coord = detCoord;
nirs.link = [link, ones(size(link,1),1)]; % append wavelength column set to 1 for CW
nirs.source.num = (1:size(srcCoord,1))';
nirs.meas.num = (1:size(detCoord,1))';

if isempty(opts.fwhm)
nirs.source.fwhm = zeros(size(srcCoord,1),1);
else
nirs.source.fwhm = opts.fwhm;
end

nirs.source.fixed = 0;
nirs.source.distributed = 0;
nirs.meas.fixed = 0;

if isscalar(opts.ri)
nirs.ri = opts.ri * ones(size(nirs.ri));
else
nirs.ri = opts.ri;
end

save_mesh(nirs, outPath);
end