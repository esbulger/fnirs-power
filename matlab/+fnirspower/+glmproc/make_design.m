function X = make_design(time, Fs, stimvec, acc_reg, extra_regs, varargin)
%MAKE_DESIGN  Build GLM design matrix consistent with the main pipeline.
%  X = glmproc.make_design(time, Fs, stimvec, acc_reg, extra_regs, ...)
%
% Positional inputs:
%   time       : T×1 (seconds). Pass [] to auto-construct from Fs.
%   Fs         : scalar Hz
%   stimvec    : T×1 binary events
%   acc_reg    : T×1 auxiliary regressor (optional, default [])
%   extra_regs : T×R additional regressors (optional, default [])
%
% Name-value pairs:
%   'hrf_seconds'    (double, default 32)
%   'normalize_cols' (logical, default true)

% -----------------------
% inputParser
% -----------------------
p = inputParser;
p.FunctionName = 'glmproc.make_design';

addRequired(p, 'time',    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addRequired(p, 'Fs',      @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addRequired(p, 'stimvec', @(x) isnumeric(x) && isvector(x));

% allow acc_reg / extra_regs to be omitted (pass [] or don’t pass them)
if nargin < 4, acc_reg = []; end
if nargin < 5, extra_regs = []; end

addOptional(p, 'acc_reg',    acc_reg,    @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addOptional(p, 'extra_regs', extra_regs, @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));

addParameter(p, 'hrf_seconds', 32, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p, 'normalize_cols', true, @(x) islogical(x) && isscalar(x));

parse(p, time, Fs, stimvec, acc_reg, extra_regs, varargin{:});
R = p.Results;

stimvec = R.stimvec(:);
T = numel(stimvec);

% Time vector
if isempty(R.time)
  time = (0:T-1)'/R.Fs;
else
  time = R.time(:);
  if numel(time) ~= T
    error('time length (%d) must match stimvec length (%d).', numel(time), T);
  end
end

% 1) trend + intercept
trend  = (1:T)';      % samples
interc = ones(T,1);

% 2) canonical HRF convolution
h = canonicalHRF(R.hrf_seconds, R.Fs);
hconv = filter(h, 1, stimvec);

% 3) assemble columns
cols = [trend, interc, hconv];

% accelerometer / aux regressor
if ~isempty(R.acc_reg)
  acc = R.acc_reg(:);
  if numel(acc) ~= T
    error('acc_reg length (%d) must match stimvec length (%d).', numel(acc), T);
  end
  acc(~isfinite(acc)) = 0; % harden

  mu = mean(acc, 'omitnan');
  sd = std(acc,  'omitnan');
  if ~isfinite(sd) || sd < eps, sd = 1; end
  accz = (acc - mu) / sd;

  cols = [cols, accz];
end

% extra regressors
if ~isempty(R.extra_regs)
  Z = R.extra_regs;
  if size(Z,1) ~= T
    error('extra_regs must have T rows (%d). Got %d.', T, size(Z,1));
  end
  Z(~isfinite(Z)) = 0;

  mu = mean(Z, 1, 'omitnan');
  sd = std(Z,  [], 1, 'omitnan');
  sd(~isfinite(sd) | sd<eps) = 1;
  Zz = (Z - mu) ./ sd;

  cols = [cols, Zz];
end

% 4) normalize columns (except intercept) to unit max
if R.normalize_cols
  scale = max(abs(cols), [], 1);
  scale(~isfinite(scale) | scale<eps) = 1;
  scale(2) = 1; % keep intercept unscaled
  cols = cols ./ scale;
end

% Final guard
cols(~isfinite(cols)) = 0;
X = cols;
end