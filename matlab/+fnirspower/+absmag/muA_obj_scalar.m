function err = muA_obj_scalar(x, mean_mu_a760, mean_mu_a850, thres_HbO, thres_HbR, ...
                              epsilon_mm_uM, J_chrom_brain, meas_HbO, meas_HbR, lambda, relative_beta_SNR)
%MUA_OBJ_SCALAR  Scalar-only wrapper for fmincon.
%  err = muA_obj_scalar(...)
%
%  Ensures the fmincon objective returns a single finite scalar.

  if nargin < 11 || isempty(relative_beta_SNR)
    relative_beta_SNR = 1.0;
  end

  err = fnirspower.absmag.muA_ls_obj_HbT_rel( ...
        x, mean_mu_a760, mean_mu_a850, thres_HbO, thres_HbR, ...
        epsilon_mm_uM, J_chrom_brain, meas_HbO, meas_HbR, lambda, relative_beta_SNR);

  % Guard: fmincon requires finite scalar
  if ~isscalar(err) || ~isfinite(err)
    err = Inf;
  end
end