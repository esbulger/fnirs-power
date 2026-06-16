function epsilon = readExtinctionCoeffs(excoefTxtPath, lambdas)
%READEXTINCTIONCOEFFS Read 2-column HbO/HbR extinction coefficients.
% excoefTxtPath: text file with columns [lambda nm, eps_HbO, eps_HbR]
% lambdas : vector of wavelengths to extract (e.g., [760 850])
% Returns epsilon as size numel(lambdas) x 2

T = importdata(excoefTxtPath);
if isstruct(T)
data = T.data;
else
data = T;
end

epsilon = nan(numel(lambdas),2);
for i = 1:numel(lambdas)
row = data(abs(data(:,1)-lambdas(i))<1e-6, 2:3);
if isempty(row)
error('Wavelength %g nm not found in %s', lambdas(i), excoefTxtPath);
end
epsilon(i,:) = row(1,:);
end
end