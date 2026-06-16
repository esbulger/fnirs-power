function [J_chrom, J_chrom_brain, pathlength_tot, pathlength_brain, J_chrom_full] = computeJChromBrain(model, epsilon_mm_uM, brain_nodes_idx)
% computeJChromBrain  Compute normalized spectral Jacobians and pathlengths
% 
%   [J_chrom, J_chrom_brain, pathlength_tot, pathlength_brain, J_chrom_full]
%     = computeJChromBrain(nirs_model, epsilon_mm_uM, brain_nodes_idx)
%   
%   Computes the pathlength-normalized spectral Jacobian (dOD/d[Hb]) for
%   a 5-layer head model, and extracts the subset corresponding to brain
%   nodes.
%
%
%   Inputs:
%     nirs_model      - struct with fields:
%                         .J_760.complex : complex Jacobian at 760 nm
%                         .J_850.complex : complex Jacobian at 850 nm
%                         .DATA760.complex: complex fluence at 760 nm
%                         .DATA850.complex: complex fluence at 850 nm
%     epsilon_mm_uM   - 2×2 extinction coefficient matrix [
%                         ε_HbO_760, ε_HbR_760;
%                         ε_HbO_850, ε_HbR_850 ] (µM^{-1}·mm^{-1})
%     brain_nodes_idx - indices of mesh nodes belonging to brain region
%
%   Outputs:
%     J_chrom         - (2N_ch × 2N_v) spectral Jacobian normalized by
%                       L-vector (total optical pathlength)
%     J_chrom_brain   - subset of J_chrom for brain nodes only (columns
%                       selected by brain_nodes_idx)
%     pathlength_tot  - total two-way pathlength vector (length 2N_ch)
%     pathlength_brain- brain-only two-way pathlength (length 2N_ch)
%     J_chrom_full    - raw (unnormalized) spectral Jacobian
%
%   Dependencies:
%     None beyond basic MATLAB functions.

%% 1. Unpack model fields
J_760 = model.J_760;       % complex Jacobian @760nm
J_850 = model.J_850;       % complex Jacobian @850nm
D760  = model.DATA760;     % complex fluence @760nm
D850  = model.DATA850;     % complex fluence @850nm
idx   = brain_nodes_idx;        % brain node indices

%% 2. Determine dimensions
n_ch       = size(J_760.complex, 1);  % number of channels
n_vertices = size(J_760.complex, 2);  % number of mesh nodes

%% 3. Compute log-intensity Jacobians (d ln I / dµ) - same as NIRFAST basic Jacobian...
phi760  = D760.complex;
I760    = abs(phi760).^2;                   % intensity @760nm
JlnI760 = 2 * real(conj(phi760)./I760 .* J_760.complex);

phi850  = D850.complex;
I850    = abs(phi850).^2;                   % intensity @850nm
JlnI850 = 2 * real(conj(phi850)./I850 .* J_850.complex);

%% 4. Assemble absorption Jacobian matrix J_mu
% Upper block for 760, lower block for 850
J_mu = [abs(JlnI760),                 zeros(n_ch, n_vertices);
        zeros(n_ch, n_vertices), abs(JlnI850)];

%% 5. Compute spectral Jacobian (unnormalized)
% Use Kronecker to apply inverse extinction matrix across channels
J_chrom_full = kron(inv(epsilon_mm_uM), eye(n_ch)) * J_mu;

%% 6. Normalize by total pathlength (L-vector)
pathlength_tot = sum(J_mu, 2) / 2;                % two-way pathlength -> one-way
J_chrom        = bsxfun(@rdivide, J_chrom_full, pathlength_tot);

%% 7. Extract brain-only columns from normalized Jacobian
% Columns for HbO then HbR at brain nodes
J_chrom_brain = J_chrom(:, [idx, idx + n_vertices]);

%% 8. Compute brain-only pathlength vector
J_mu_brain      = [abs(JlnI760(:, idx)), zeros(n_ch, numel(idx));
                   zeros(n_ch, numel(idx)), abs(JlnI850(:, idx))];
pathlength_brain= sum(J_mu_brain, 2) / 2;

end
