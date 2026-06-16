function E = epochs_from_arrays(all_hbo, all_hbr, Fs, baseline_sec, block_sec, all_aux)
%EPOCHS_FROM_ARRAYS Build epoch struct E from epoched arrays in memory.
% E = nirsproc.epochs_from_arrays(all_hbo, all_hbr, Fs, baseline_sec, block_sec, all_aux)
%
% Inputs
% all_hbo : nBlocks × T × nCh
% all_hbr : nBlocks × T × nCh
% Fs : sampling rate
% baseline_sec, block_sec : stored in E for downstream design building
% all_aux (optional): nBlocks × T real aux epochs. If omitted, proxy used.
%
% Output
% E : struct with fields .hbo, .hbr, .aux (proxy or real), .aux_real (optional), .Fs, .T


E.hbo = all_hbo; E.hbr = all_hbr; E.Fs = Fs; E.T = size(all_hbo,2);

if nargin >= 6 && ~isempty(all_aux)
E.aux = all_aux; E.aux_real = all_aux;
else
E.aux = squeeze(mean(abs(all_hbo), 3));
end
E.baseline_sec = baseline_sec; E.block_sec = block_sec;
end