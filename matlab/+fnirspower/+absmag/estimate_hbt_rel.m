function result = estimate_hbt_rel(inputs)
%ESTIMATE_HBT_REL  Estimate ΔHbO/ΔHbR using spectral Jacobian and masks.
%  result = ESTIMATE_HBT_REL(inputs)
%
%  Inputs (struct)
%    .mean_mu_a760, .mean_mu_a850 : brain-only µa templates (nBrain×1)
%    .J_chrom_brain               : spectral Jacobian (2*nCh × 2*nBrain)
%    .meas_HbO, .meas_HbR         : measured betas (nCh×1)
%    .thres_HbO, .thres_HbR       : logical masks (1×nCh)
%    .epsilon_mm_uM               : 2×2 epsilon matrix
%    .relative_beta_SNR           : scalar
%    .lambda                      : scalar (regularization weight)
%
%  Output struct result
%    .x_opt (ΔHbO, ΔHbR), .err_tot, .err_O, .err_R, .err_rat
%    .dmu_a_opt (2×1), .HbO_est, .HbR_est

args = inputs;

objFun = @(x) fnirspower.absmag.muA_ls_obj_HbT_rel(x, args.mean_mu_a760, args.mean_mu_a850, ...
                          args.thres_HbO, args.thres_HbR, args.epsilon_mm_uM, ...
                          args.J_chrom_brain, args.meas_HbO, args.meas_HbR, ...
                           args.lambda, args.relative_beta_SNR);

opts = optimoptions('fmincon','Display','none','OptimalityTolerance',1e-8);
[x_opt, ~] = fmincon(@(x) objFun(x), [5.70; -1.91], [],[],[],[], [0;-Inf], [Inf;0], [], opts);
[err_tot, err_O, err_R, err_rat, dHb_opt, dmu_a_opt] = objFun(x_opt);

% Build estimated signals
m760 = args.mean_mu_a760 / max(abs(args.mean_mu_a760));
m850 = args.mean_mu_a850 / max(abs(args.mean_mu_a850));
vec_norm = [m760 * dmu_a_opt(1); m850 * dmu_a_opt(2)];
Delta_Hb = args.J_chrom_brain * vec_norm;
nCh = numel(args.meas_HbO);
HbO_est = Delta_Hb(1:nCh);
HbR_est = Delta_Hb(nCh+1:end);
HbO_est(args.thres_HbO) = 0; 
HbR_est(args.thres_HbR) = 0;

result = struct('x_opt', x_opt, 'err_tot',err_tot, 'err_O',err_O, 'err_R',err_R, 'err_rat',err_rat, ...
                'dmu_a_opt', dmu_a_opt, 'HbO_est', HbO_est, 'HbR_est', HbR_est);
end
