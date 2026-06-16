function [err, err_O, err_R, err_rat, dHb_opt, dmu_a_opt] = muA_ls_obj_HbT_rel( ...
        x, mean_mu_a760, mean_mu_a850, thres_HbO, thres_HbR, ...
        epsilon_mm_uM, J_chrom_brain, meas_HbO, meas_HbR, lambda, relative_beta_SNR)
%MUA_LS_OBJ_HBT_REL  HbO/HbR LS fit + 3:1 ratio regularization (+ optional SNR weighting).
%
%  x = [ΔHbO; ΔHbR] (µM), with constraints ΔHbO≥0, ΔHbR≤0 (typically enforced in fmincon bounds)
%  lambda : nonnegative weight on the ratio penalty (toward r0=3)
%  relative_beta_SNR : scalar SNR(HbO)/SNR(HbR). HbR term is weighted by 1/relative_beta_SNR.
%
%  Returns:
%    err      : total scalar objective
%    err_O    : HbO misfit (∑(meas_O - pred_O)^2)
%    err_R    : HbR misfit (∑(meas_R - pred_R)^2)
%    err_rat  : ratio penalty ( (ΔHbO/|ΔHbR| - r0)^2 )
%    dHb_opt  : [ΔHbO; ΔHbR]
%    dmu_a_opt: [Δμa@760; Δμa@850]
%

  % ---- Optional argument ----
  if nargin < 12 || isempty(relative_beta_SNR)
    relative_beta_SNR = 1.0;
  end

  % Encouraged HbO:|HbR| ratio - to prevent unrealistic result
  r0 = 3;

  % Store Hb
  dHb_opt = x(:);                    % [HbO; HbR]

  % Convert to absorption changes at 760/850
  dmu_a_opt = epsilon_mm_uM * dHb_opt; % [dμa760; dμa850]

  % Normalize voxelwise patterns to unit peak (avoid scale effects)
  m760 = mean_mu_a760 / max(eps, max(abs(mean_mu_a760)));
  m850 = mean_mu_a850 / max(eps, max(abs(mean_mu_a850)));

  % Build stacked μa vector matching J_chrom_brain column ordering (760 brain; 850 brain)
  vec_muA_norm = [m760 * dmu_a_opt(1); m850 * dmu_a_opt(2)];  % (2*nBrain)×1

  % Forward to channels in chromophore space
  Delta_Hb = J_chrom_brain * vec_muA_norm;   % (2*nCh)×1
  nCh = numel(meas_HbO);

  HbO_pred = Delta_Hb(1:nCh)';
  HbR_pred = Delta_Hb(nCh+1:end)';

  % Apply masks (zero both pred & meas where thresholded)
  if ~isempty(thres_HbO)
    HbO_pred(thres_HbO) = 0;  meas_HbO(thres_HbO) = 0;
  end
  if ~isempty(thres_HbR)
    HbR_pred(thres_HbR) = 0;  meas_HbR(thres_HbR) = 0;
  end

  % Data misfit (scalar)
  err_O = sum((meas_HbO - HbO_pred).^2);
  err_R = sum((meas_HbR - HbR_pred).^2);

  % Ratio penalty toward r0 = 3
  rat     = x(1) / max(eps, abs(x(2)));
  err_rat = (rat - r0)^2;

  % Composite objective (HbR down-weighted by relative SNR)
  wR  = 1 / max(eps, relative_beta_SNR);
  err = err_O + wR * err_R + lambda * err_rat;
end