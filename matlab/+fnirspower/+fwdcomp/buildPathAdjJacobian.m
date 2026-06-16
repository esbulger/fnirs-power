function SJ = buildPathAdjJacobian(nirs_model, epsilon_mm_uM, brain_nodes_idx)
%BUILD_SPECTRAL_JACOBIAN  Compute spectral Jacobians and pathlengths.
%  SJ = BUILD_SPECTRAL_JACOBIAN(nirs_model, epsilon_mm_uM, brain_nodes_idx)
%
%  Output struct SJ
%    .J_chrom         : (2*nCh)×(2*nV) chromophore Jacobian
%    .J_chrom_brain   : (2*nCh)×(2*nBrain) brain-only columns
%    .pathlength_tot  : (2*nCh)×1 total pathlengths (760 then 850)
%    .pathlength_brain: (2*nCh)×1 brain-only pathlengths
%
% This version is an effective pathlength modified version, intended for
% prediction of real-valued channel-level measurements according to SD
% distance, as is implicitly computed in NIRFAST
%

[J_chrom, J_chrom_brain, pathlength_tot, pathlength_brain] = ...
  fnirspower.fwdcomp.computeJChromBrain(nirs_model, epsilon_mm_uM, brain_nodes_idx);

SJ = struct('J_chrom', J_chrom, 'J_chrom_brain', J_chrom_brain, ...
            'pathlength_tot', pathlength_tot, 'pathlength_brain', pathlength_brain);
end
