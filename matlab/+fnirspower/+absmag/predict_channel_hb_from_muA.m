function [HbO_pred, HbR_pred, Hb_all] = predict_channel_hb_from_muA(J_chrom_brain, mu_a760, mu_a850)
%PREDICT_CHANNEL_HB_FROM_MUA  Apply spectral J to brain µa vectors → channels.
%  [HbO_pred, HbR_pred, Hb_all] = PREDICT_CHANNEL_HB_FROM_MUA(Jc_brain, mu760, mu850)

mu_vec = [mu_a760(:); mu_a850(:)];
Hb_all = J_chrom_brain * mu_vec;
nCh = numel(Hb_all)/2;
HbO_pred = Hb_all(1:nCh).';
HbR_pred = Hb_all(nCh+1:end).';
end
