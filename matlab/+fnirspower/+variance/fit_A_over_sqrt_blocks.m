function fit_out = fit_A_over_sqrt_blocks(iter_blocks, ydata)
%FIT_A_OVER_SQRT_BLOCKS Fit y = A/sqrt(B) to block-count data.
%
%  fit_out = fnirspower.variance.fit_A_over_sqrt_blocks(iter_blocks, ydata)
%
%  Inputs
%  ------
%  iter_blocks : numeric vector
%      Block counts B.
%
%  ydata : numeric vector
%      Observed values at each block count.
%
%  Returns
%  -------
%  fit_out : struct with fields
%      .A       - fitted scalar coefficient
%      .y_fit   - fitted values at iter_blocks
%      .R2      - coefficient of determination
%      .xdata   - block counts used for fit
%      .ydata   - observed values used for fit

arguments
    iter_blocks (1,:) double
    ydata (1,:) double
end

xdata = iter_blocks(:);
ydata = ydata(:);

valid = isfinite(xdata) & isfinite(ydata) & xdata > 0;
xdata = xdata(valid);
ydata = ydata(valid);

if isempty(xdata)
    error('No valid data points available for A/sqrt(B) fit.');
end

model_fun = @(A,x) A ./ sqrt(x);
A0 = 1;
A_fit = lsqcurvefit(model_fun, A0, xdata, ydata);

y_fit = model_fun(A_fit, xdata);
SSres = sum((ydata - y_fit).^2);
SStot = sum((ydata - mean(ydata)).^2);
R2 = 1 - SSres / max(eps, SStot);

fit_out = struct();
fit_out.A = A_fit;
fit_out.y_fit = y_fit(:).';
fit_out.R2 = R2;
fit_out.xdata = xdata(:).';
fit_out.ydata = ydata(:).';
end