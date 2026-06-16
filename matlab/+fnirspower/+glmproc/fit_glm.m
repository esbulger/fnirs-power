function betas = fit_glm(Y, X, Fs, opts)
%GLMPROC.FIT_GLM  AR-IRLS GLM fit for multichannel time series.
%  betas = GLMPROC.FIT_GLM(Y, X, Fs, opts)
%
%  Inputs
%    Y    : T×M matrix (time × channels)
%    X    : T×p design matrix
%    Fs   : sampling rate (Hz)
%    opts : struct with field
%             .Pmax : AR order for prewhitening (default ceil(4*Fs))
%
%  Output
%    betas : p×M matrix of GLM coefficients; use betas(3,:) for task reg.
%
%  Requires: glm_ar_irls from Huppert NIRS toolbox.

arguments
  Y double
  X double
  Fs (1,1) double
  opts.Pmax double = []
end

if isempty(opts.Pmax), Pmax = ceil(4*Fs); else, Pmax = opts.Pmax; end
res = glm_ar_irls(Y, X, Pmax);
betas = res; % return all for flexibility; caller can select row 3.
end
