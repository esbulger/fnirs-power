function muA_file = reconstruct_mu_a(beta_od_file, nirsmodel_file, mesh_mat, save_dir, tik_alpha)
%FNIRSPOWER.RECON.RECONSTRUCT_MU_A
% Reconstruct subject-level absorption-change maps at 760 and 850 nm.
%
% This function reconstructs voxelwise absorption changes from channel-level
% optical-density beta estimates using Tikhonov-regularized inversion. The
% full-volume reconstructions are restricted to grey-matter nodes.
%
% Usage
% -----
% muA_file = fnirspower.recon.reconstruct_mu_a( ...
%     beta_od_file, nirsmodel_file, mesh_mat, save_dir, tik_alpha);
%
% Inputs
% ------
% beta_od_file :
%     Path to a MAT file containing:
%
%       all_beta_a760
%           Optical-density beta estimates at 760 nm
%           [nSubjects x nChannels].
%
%       all_beta_a850
%           Optical-density beta estimates at 850 nm
%           [nSubjects x nChannels].
%
% nirsmodel_file :
%     Path to a forward-model MAT file containing J_760 and J_850.
%
% mesh_mat :
%     Path to a mesh MAT file containing genmesh.ele. Tissue region 2 is
%     treated as grey matter.
%
% save_dir :
%     Directory in which to save the reconstructed absorption maps. If
%     empty, the output is saved beside beta_od_file.
%
% tik_alpha :
%     Tikhonov regularization weight passed to
%     fnirspower.recon.invert_woodbury.
%
% Output
% ------
% muA_file :
%     Full path to the saved MAT file. The file contains:
%
%       all_mu_a760
%           Reconstructed absorption changes at 760 nm
%           [nSubjects x nGreyMatterNodes].
%
%       all_mu_a850
%           Reconstructed absorption changes at 850 nm
%           [nSubjects x nGreyMatterNodes].
%
%       brain_nodes_idx
%           Indices of grey-matter nodes in the full mesh.
%
% Notes
% -----
% Channels with nonfinite beta estimates are excluded separately for each
% wavelength during reconstruction.
%
% See also:
%   fnirspower.recon.invert_woodbury

% Load ΔOD betas
B = load(beta_od_file);
if ~isfield(B,'all_beta_a760') || ~isfield(B,'all_beta_a850')
  error('ΔOD MAT must contain all_beta_a760 and all_beta_a850 (nSubj×nCh).');
end
all_beta_a760 = B.all_beta_a760;        % nSubj×nCh
all_beta_a850 = B.all_beta_a850;        % nSubj×nCh
nSubj = size(all_beta_a760,1);

% Brain indices from mesh
M = load(mesh_mat,'genmesh'); e = M.genmesh.ele;
brain_tet = e(e(:,5)==2,:);
brain_nodes_idx = unique([brain_tet(:,1); brain_tet(:,2); brain_tet(:,3); brain_tet(:,4)]);

% NIRFAST Jacobians
NM = load(nirsmodel_file,'J_760','J_850');
J_760 = abs(NM.J_760.complete);
J_850 = abs(NM.J_850.complete);

alpha = tik_alpha;  % Tikhonov weight

all_mu_a760 = NaN(nSubj, numel(brain_nodes_idx));
all_mu_a850 = NaN(nSubj, numel(brain_nodes_idx));

for i = 1:nSubj
  % Determine bad channels per wavelength (non-finite betas)
  y760 = all_beta_a760(i,:)'; bad760 = find(~isfinite(y760));
  y850 = all_beta_a850(i,:)'; bad850 = find(~isfinite(y850));

  % Reconstruct full-volume map, then restrict to brain nodes
  x760 = fnirspower.recon.invert_woodbury(J_760, y760, bad760, alpha);
  x850 = fnirspower.recon.invert_woodbury(J_850, y850, bad850, alpha);

  all_mu_a760(i,:) = x760(brain_nodes_idx)';
  all_mu_a850(i,:) = x850(brain_nodes_idx)';
end

% Decide save location
if ~isempty(save_dir)
  if exist(save_dir,'dir')~=7, mkdir(save_dir); end
  out_dir = save_dir;
else
  out_dir = fileparts(beta_od_file);
end

% Save recon MAT in the standard format expected downstream
date_str = datestr(datetime('now'),'yyyy-mm-dd');
muA_file = fullfile(out_dir, [date_str '_mu_a_recon.mat']);
save(muA_file,'all_mu_a760','all_mu_a850','brain_nodes_idx','-v7.3');
fprintf('Built μa recon (recon.invert_woodbury) → %s\n', muA_file);
end