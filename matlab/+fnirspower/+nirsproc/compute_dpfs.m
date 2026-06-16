function [DPF_760, DPF_850, info] = compute_dpfs(model, epsilon_mm_uM, brain_nodes_idx, wl_idx, dist_ch, nCh)
%COMPUTE_DPFS  Channel DPFs from Jacobian-derived total pathlengths.
%  [DPF_760, DPF_850, info] = nirsproc.compute_dpfs(model, epsilon_mm_uM, brain_nodes_idx, wl_idx, dist_ch, nCh)

arguments
  model struct
  epsilon_mm_uM double
  brain_nodes_idx (:,1) double
  wl_idx (:,1) double
  dist_ch (:,1) double
  nCh (1,1) double
end

% Get total pathlengths per measurement from chromophore Jacobian helper
[~, ~, pathlength_tot] = fnirspower.fwdcomp.computeJChromBrain(model, epsilon_mm_uM, brain_nodes_idx);
M = numel(pathlength_tot);

info = struct('M',M,'nCh',nCh);

% Case A: vector is exactly [760(1..nCh); 850(1..nCh)]
if M == 2*nCh
  idx760 = (1:nCh).';
  idx850 = (nCh+1:2*nCh).';

% Case B: vector matches the full measurement list — pick by wl_idx
elseif M == numel(wl_idx)
  idx760 = find(wl_idx==1);
  idx850 = find(wl_idx==2);
  % Be robust to overcomplete ordering
  if numel(idx760) < nCh || numel(idx850) < nCh
    error('compute_dpfs:idx', 'Found fewer than nCh measurements per wavelength (760=%d, 850=%d, nCh=%d).', numel(idx760), numel(idx850), nCh);
  end
  idx760 = idx760(1:nCh);
  idx850 = idx850(1:nCh);
else
  error('compute_dpfs:size', 'pathlength_tot length (%d) not equal to 2*nCh (%d) or numel(wl_idx) (%d).', M, 2*nCh, numel(wl_idx));
end

% Each is nCh×1
DPF_760 = pathlength_tot(idx760) ./ dist_ch(:);  % pathlength_tot already adjusted for 2-way pathlength
DPF_850 = pathlength_tot(idx850) ./ dist_ch(:);

info.idx760 = idx760; info.idx850 = idx850;
end
