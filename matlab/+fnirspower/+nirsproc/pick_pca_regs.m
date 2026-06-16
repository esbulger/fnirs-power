function [pca_regs_Hb, pca_regs_OD, diag_info] = pick_pca_regs(HBO, HBD, dOD_hp, conv, n_pca, r2max)
%NIRSPROC.PICK_PCA_REGS  Select task-uncorrelated PCs as nuisance regressors.
%  [pca_regs_Hb, pca_regs_OD, diag_info] = NIRSPROC.PICK_PCA_REGS(...)
%
%  Inputs
%    HBO, HBD : T×nCh HbO/HbR time series
%    dOD_hp   : T×(2*nCh) ΔOD time series (high-passed)
%    conv     : T×1 task regressor (stim * HRF)
%    n_pca    : number of top PCs to consider (default 1)
%    r2max    : percent R² threshold with conv to allow inclusion (default 2.5)
%
%  Outputs
%    pca_regs_Hb : T×K PCs from Hb data that pass R² criterion
%    pca_regs_OD : T×K PCs from ΔOD that pass R² criterion
%    diag_info   : struct array with fields {which, pc_index, r, R2, var_expl}
%
%  Author: Eli Bulger (refactor)
%  Version: 0.1.0

arguments
  HBO double
  HBD double
  dOD_hp double
  conv (:,1) double
  n_pca (1,1) double = 1
  r2max (1,1) double = 2.5
end

Hb = [HBO, HBD];
[pca_regs_Hb, diag_Hb] = pick_from_matrix(Hb, conv, n_pca, r2max, 'Hb');
[pca_regs_OD, diag_OD] = pick_from_matrix(dOD_hp, conv, n_pca, r2max, 'OD');

diag_info = [diag_Hb, diag_OD];
end

function [regs, diag_info] = pick_from_matrix(M, conv, n_pca, r2max, which)
  [~,scores,~,~,expl] = pca(M);
  regs = []; diag_info = struct('which',{},'pc_index',{},'r',{},'R2',{},'var_expl',{});
  for jj=1:n_pca
    pc = scores(:,jj); pc = pc ./ max(1,abs(max(pc)));
    r = corr(pc, conv); R2 = 100*r^2; ve = expl(jj);
    if R2 <= r2max
      regs = [regs, pc];
    end
    diag_info(end+1) = struct('which',which,'pc_index',jj,'r',r,'R2',R2,'var_expl',ve);
  end
end
