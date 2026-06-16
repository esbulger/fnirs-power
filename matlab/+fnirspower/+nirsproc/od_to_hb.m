function [HBO, HBD] = od_to_hb(dOD_hp, nCh, DPF760, DPF850, dist_ch, coeffs)
%NIRSPROC.OD_TO_HB  Convert ΔOD (two wavelengths) to HbO/HbR time series.
%  [HBO, HBD] = NIRSPROC.OD_TO_HB(dOD_hp, nCh, DPF760, DPF850, dist_ch, coeffs)
%
%  Inputs
%    dOD_hp  : T×(2*nCh) matrix, high-passed ΔOD (760 then 850 columns)
%    nCh     : number of channels (per wavelength)
%    DPF760  : nCh×1 differential pathlength factor for 760 nm
%    DPF850  : nCh×1 differential pathlength factor for 850 nm
%    dist_ch : nCh×1 source–detector distance (mm)
%    coeffs  : struct with extinction coefficients (mm^-1/mM)
%              .hbo760 .hbo850 .hbr760 .hbr850
%
%  Outputs
%    HBO, HBD : T×nCh in µM
%
%  Notes
%    Uses per-channel 2×2 forward matrix and pinv. Multiplied by 1000 to
%    convert mM to µM, matching upstream convention.
%
%  Author: Eli Bulger (refactor)
%  Version: 0.1.0

arguments
  dOD_hp double
  nCh (1,1) double
  DPF760 (:,1) double
  DPF850 (:,1) double
  dist_ch (:,1) double
  coeffs struct
end

HBO = NaN(size(dOD_hp,1), nCh); HBD = HBO;
for j=1:nCh
  F = [DPF760(j)*coeffs.hbo760*dist_ch(j), DPF760(j)*coeffs.hbr760*dist_ch(j);
       DPF850(j)*coeffs.hbo850*dist_ch(j), DPF850(j)*coeffs.hbr850*dist_ch(j)];
  dod_760 = dOD_hp(:, j); dod_850 = dOD_hp(:, j+nCh);
  hb = 1000 * [dod_760, dod_850] * pinv(F)'; % µM
  HBO(:,j) = hb(:,1); HBD(:,j) = hb(:,2);
end
end
