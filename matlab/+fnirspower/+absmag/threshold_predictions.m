function [HbO_thr, HbR_thr, thres_HbO, thres_HbR] = threshold_predictions(meas_HbO, meas_HbR, all_HbO_pred, all_HbR_pred, alpha)
%THRESHOLD_PREDICTIONS  t-test across subjects + align with measured masks.
%  [HbO_thr, HbR_thr, thres_HbO, thres_HbR] = THRESHOLD_PREDICTIONS(...)

if nargin<5 || isempty(alpha), alpha = 0.05; end

[~,p_O] = ttest(all_HbO_pred, 0, 'Alpha', alpha);
[~,p_R] = ttest(all_HbR_pred, 0, 'Alpha', alpha);

thres_meas_HbO = (meas_HbO == 0);
thres_meas_HbR = (meas_HbR == 0);
thres_fwd_HbO  = (p_O > alpha);
thres_fwd_HbR  = (p_R > alpha);

thres_HbO = thres_meas_HbO(:).' | thres_fwd_HbO(:).';
thres_HbR = thres_meas_HbR(:).' | thres_fwd_HbR(:).';

HbO_thr = mean(all_HbO_pred,1);
HbR_thr = mean(all_HbR_pred,1);
HbO_thr(thres_HbO) = 0;
HbR_thr(thres_HbR) = 0;
end
