function nirs = setOpticalProps(nirs, props)
%SETOPTICALPROPS Fill mua/mus by region for a given wavelength.
% props is a struct with region fields 1..5 each containing .mua and .mus (mm^-1)
% Example: props(760).r(5).mua = 0.0170; props(760).r(5).mus = 0.74; etc.

regions = unique(nirs.region(:));
for rid = regions.'
if ~isfield(props, 'r') || numel(props.r) < rid
error('Missing optical properties for region %d', rid);
end
nirs.mua(nirs.region==rid) = props.r(rid).mua;
nirs.mus(nirs.region==rid) = props.r(rid).mus;
end
end