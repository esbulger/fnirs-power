function dist_ch = sd_distances(raw, link)
%NIRSPROC.SD_DISTANCES  Compute source–detector distances (mm) per channel.
%  dist_ch = SD_DISTANCES(raw, link)

nCh = size(link,1);
dist_ch = zeros(nCh,1);
for j=1:nCh
  dist_ch(j) = norm(raw.nirs.probe.sourcePos3D(link(j,1),:) - ...
                   raw.nirs.probe.detectorPos3D(link(j,2),:));
end
end
