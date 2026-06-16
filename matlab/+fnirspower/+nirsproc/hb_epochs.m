function E = hb_epochs(S, baseline_sec, block_sec, aux_timeseries)
%HB_EPOCHS Epoch HbO/HbR (and aux) from a preprocessed subject.
% E = nirsproc.hb_epochs(S, baseline_sec, block_sec, aux_timeseries)
%
% Inputs
% S : struct from nirsproc.preprocess (one subject)
% baseline_sec : seconds before onset to include in each epoch
% block_sec : seconds after onset (task duration) to include
% aux_timeseries (optional) : T_full×1 accel/aux signal to epoch.
% If omitted, uses proxy = mean |dOD_hp| across channels.
%
% Outputs (struct E)
% .hbo : nBlocks × T × nCh
% .hbr : nBlocks × T × nCh
% .aux : nBlocks × T (will equal .aux_real if provided)
% .aux_real : nBlocks × T (present only if aux_timeseries provided)
% .Fs, .T, .idx_onsets, .baseline_sec, .block_sec


arguments
S struct
baseline_sec (1,1) double
block_sec (1,1) double
aux_timeseries double = []
end


Fs = S.Fs; on = find(S.stimvec>0.5);
if isempty(on), error('HB_EPOCHS: no stimulus onsets found'); end
on = on([true; diff(on)>round(4*Fs)]); % ensure spacing


pre = round(baseline_sec*Fs);
post = round((block_sec+15)*Fs); % +15 s post like legacy
T = pre+post;


nBlocks = numel(on);
nCh = S.nCh;
E.hbo = NaN(nBlocks, T, nCh);
E.hbr = NaN(nBlocks, T, nCh);
E.aux = NaN(nBlocks, T);


% default aux proxy
aux_proxy = mean(abs(S.dOD_hp), 2);
use_real_aux = ~isempty(aux_timeseries) && numel(aux_timeseries) == numel(S.time);


for b = 1:nBlocks
idx = on(b);
a = idx-pre+1; z = idx+post; % inclusive
if a<1 || z>numel(S.time), continue; end
E.hbo(b,:,:) = S.HBO(a:z, :);
E.hbr(b,:,:) = S.HBD(a:z, :);
if use_real_aux
E.aux(b,:) = aux_timeseries(a:z,1).';
else
E.aux(b,:) = aux_proxy(a:z,1).';
end
end


if use_real_aux
E.aux_real = E.aux;
end


E.Fs = Fs; E.T = T; E.idx_onsets = on;
E.baseline_sec = baseline_sec; E.block_sec = block_sec;
end