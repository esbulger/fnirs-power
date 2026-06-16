function X = build_concat_design(Fs, block_T, nBlocks, baseline_sec, block_sec, aux, varargin)
%BUILD_CONCAT_DESIGN  Design matrix for concatenated epochs via glmproc.make_design.
%  X = variance.build_concat_design(Fs, block_T, nBlocks, baseline_sec, block_sec, aux, ...)
%
%  Positional inputs:
%    Fs           : scalar Hz
%    block_T      : samples per block (epoch length in samples)
%    nBlocks      : number of blocks concatenated
%    baseline_sec : baseline duration (s) inside each epoch
%    block_sec    : task duration (s) inside each epoch
%    aux          : (Ttot×1) concatenated auxiliary regressor (optional; pass [] to omit)
%
%  Name-value pairs:
%    'hrf_seconds'    : default 32
%    'normalize_cols' : default true
%
%  Notes:
%    - Builds a tiled stimvec for the concatenated series and calls
%      fnirspower.glmproc.make_design for consistency.
%    - Sanitizes aux to avoid NaN/Inf propagating into design matrix.

% -----------------------
% inputParser
% -----------------------
p = inputParser;
p.FunctionName = 'variance.build_concat_design';

addRequired(p, 'Fs',           @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addRequired(p, 'block_T',      @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addRequired(p, 'nBlocks',      @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addRequired(p, 'baseline_sec', @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0);
addRequired(p, 'block_sec',    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);

addOptional(p, 'aux', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));

addParameter(p, 'hrf_seconds', 32, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p, 'normalize_cols', true, @(x) islogical(x) && isscalar(x));

parse(p, Fs, block_T, nBlocks, baseline_sec, block_sec, aux, varargin{:});
opts = p.Results;

% -----------------------
% Build stimvec
% -----------------------
Ttot = round(block_T) * round(nBlocks);

stim_block = zeros(round(block_T), 1);
i0 = max(1, round(opts.baseline_sec * opts.Fs));

on  = round(baseline_sec * Fs);
off = round((baseline_sec + block_sec) * Fs);

on  = max(1, min(round(block_T), on));
off = max(1, min(round(block_T), off));

stim_block(on:off) = 1;
stimvec = repmat(stim_block, round(nBlocks), 1);

time = (0:Ttot-1)'/Fs;

% -----------------------
% Aux handling (robust)
% -----------------------
acc_reg = [];
if ~isempty(opts.aux)
  acc_reg = opts.aux(:);

  % enforce length match
  if numel(acc_reg) ~= Ttot
    error('aux length (%d) does not match expected concatenated length Ttot (%d).', numel(acc_reg), Ttot);
  end

  % Replace NaN/Inf safely
  bad = ~isfinite(acc_reg);
  if any(bad)
    % fill with median of finite values, else zeros
    finite_vals = acc_reg(isfinite(acc_reg));
    if isempty(finite_vals)
      acc_reg(:) = 0;
    else
      acc_reg(bad) = median(finite_vals);
    end
  end
end

% Call canonical design builder
X = fnirspower.glmproc.make_design( ...
      time, Fs, stimvec, acc_reg, [], ...
      'hrf_seconds', opts.hrf_seconds, ...
      'normalize_cols', opts.normalize_cols);

% Final guard: no NaN/Inf in X
if any(~isfinite(X(:)))
  X(~isfinite(X)) = 0;
end
end