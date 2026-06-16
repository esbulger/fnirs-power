function [dOD_760, dOD_850] = split_wavelengths(dOD_hp, wl_idx)
%NIRSPROC.SPLIT_WAVELENGTHS  Split ΔOD columns into 760 and 850 blocks.
%  [dOD_760, dOD_850] = SPLIT_WAVELENGTHS(dOD_hp, wl_idx)

dOD_760 = dOD_hp(:, wl_idx==1);
dOD_850 = dOD_hp(:, wl_idx==2);
end
