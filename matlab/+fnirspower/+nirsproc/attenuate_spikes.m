function X = attenuate_spikes(X, win, min_peak_dist, min_peak_height)
%NIRSPROC.ATTENUATE_SPIKES  Heuristic spike squashing via derivative energy.
%  X = ATTENUATE_SPIKES(X, win, min_peak_dist, min_peak_height)

if nargin<2, win=30; end
if nargin<3, min_peak_dist=30; end
if nargin<4, min_peak_height=0.02; end

warning('off','signal:findpeaks:largeMinPeakHeight');
for j=1:size(X,2)
  tmp = X(:,j);
  sqdiff = movsum(diff(tmp).^2, [win-1 0]); sqdiff = sqdiff(max(1,win):end);
  [~,loc] = findpeaks(sqdiff,'MinPeakDistance',min_peak_dist,'MinPeakHeight',min_peak_height);
  for ii=1:length(loc)
    a = loc(ii); b = min(a+win, numel(tmp));
    tmp(a:b) = tmp(a) + (0:(b-a))' .* (tmp(b)-tmp(a)) / max(1,(b-a));
  end
  X(:,j) = tmp;
end
end
