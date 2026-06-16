function stimvec = build_stimvec(raw, Fs, T)
%NIRSPROC.BUILD_STIMVEC  Build a T×1 impulse vector from stim onsets.
%  stimvec = BUILD_STIMVEC(raw, Fs, T)

stimvec = zeros(T,1);
for k=1:numel(raw.nirs.stim)
  if isempty(raw.nirs.stim(k).data), continue; end
  idx = round(raw.nirs.stim(k).data(:,1)*Fs);
  idx = idx(idx>=1 & idx<=T);
  stimvec(idx) = 1;
end
end
