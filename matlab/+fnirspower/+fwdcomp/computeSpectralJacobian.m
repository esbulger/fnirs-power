function Jspec = computeSpectralJacobian(Jcell, epsilon)
%COMPUTESPECTRALJACOBIAN Convert dOD/dmua to dOD/d[HbO,HbR].
% Jcell : cell array of Jacobian structs per wavelength with field .complete
% epsilon : [L x 2] extinction coefficients matching Jcell order
% Returns : [ (#meas*L) x 2N ] stacked spectral Jacobian over nodes

L = numel(Jcell);
if size(epsilon,1) ~= L
error('epsilon rows (%d) must equal number of wavelengths (%d).', size(epsilon,1), L);
end

rows = cell(L,1);
for i = 1:L
Ji = Jcell{i}.complete; % dOD/dmua
rows{i} = [Ji*epsilon(i,1), Ji*epsilon(i,2)];
end
Jspec = vertcat(rows{:});
end