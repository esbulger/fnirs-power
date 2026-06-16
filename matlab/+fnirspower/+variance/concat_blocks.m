function Y = concat_blocks(E_blocks, baseline_samples, taper_samples)
%CONCAT_BLOCKS Concatenate epochs with default stitched-baseline alignment.
% Y = variance.concat_blocks(E_blocks, baseline_samples, taper_samples)
%
% Default behavior
% ----------------
% Prepare a continuous time series from multiple block epochs by
% (1) aligning the start-baseline of each block to the end of the
%     preceding block, per channel,
% (2) applying a low-order polynomial detrend to the concatenated series
%     to reduce slow drift introduced by stitching
%
% Legacy behavior
% ---------------
% If taper_samples is passed as a numeric scalar, the function uses the
% original behavior:
% (1) per-block baseline subtraction over the first baseline_samples
%     points, (2) optional edge taper, and (3) concatenation.
%
% Optional struct behavior
% ------------------------
% taper_samples may also be a struct with fields such as:
%   opts.mode              = 'baseline_stitch'   % default
%                         or 'baseline_taper'
%   opts.taper_samples     = 0
%   opts.poly_order        = 1
%   opts.end_match_samples = []
%
% Inputs
% ------
% E_blocks : B x T x nCh (or B x T x 1 for aux)
% baseline_samples : integer > 0; number of baseline samples
%
% taper_samples :
%   - omitted or []: default stitched-baseline mode
%   - numeric scalar: legacy baseline-subtract + taper mode
%   - struct: explicit options
%
% Output
% ------
% Y : (B*T) x nCh concatenated matrix (time x channels)

arguments
    E_blocks
    baseline_samples (1,1) double {mustBePositive}
    taper_samples = []
end

B   = size(E_blocks,1);   % number of blocks
T   = size(E_blocks,2);   % samples per block
nCh = size(E_blocks,3);   % channels
X   = E_blocks;

L0 = min(T, round(baseline_samples));

% -------------------------------------------------------------------------
% Parse options while preserving backward compatibility
% -------------------------------------------------------------------------
opts = struct();
opts.mode              = 'baseline_stitch';
opts.taper_samples     = 0;
opts.poly_order        = 1;
opts.end_match_samples = [];

if nargin < 3 || isempty(taper_samples)
    % keep defaults
elseif isnumeric(taper_samples) && isscalar(taper_samples)
    % backward-compatible legacy behavior
    opts.mode          = 'baseline_taper';
    opts.taper_samples = taper_samples;
elseif isstruct(taper_samples)
    user_opts = taper_samples;
    fn = fieldnames(user_opts);
    for i = 1:numel(fn)
        opts.(fn{i}) = user_opts.(fn{i});
    end
    if ~isfield(user_opts, 'taper_samples') || isempty(user_opts.taper_samples)
        opts.taper_samples = 0;
    end
else
    error('taper_samples must be empty, a numeric scalar, or an options struct.');
end

if isempty(opts.end_match_samples)
    Lend = L0;
else
    Lend = min(T, round(opts.end_match_samples));
end

% -------------------------------------------------------------------------
% Mode 1: Per-block baseline subtraction + optional taper
% -------------------------------------------------------------------------
if strcmpi(opts.mode, 'baseline_taper')

    % === (1) Per-block baseline alignment (per channel) ===
    mu0 = squeeze(mean(X(:,1:L0,:), 2, 'omitnan')); % B x nCh
    for b = 1:B
        X(b,:,:) = X(b,:,:) - reshape(mu0(b,:), [1 1 nCh]);
    end

    % === (2) Optional cosine taper at block starts/ends ===
    if opts.taper_samples > 0
        L = min(round(opts.taper_samples), floor(T/4)); % guard over-tapering
        if L > 0
            w  = hann(2*L);
            wL = w(1:L);
            wR = w(L+1:end);
            for b = 1:B
                X(b,1:L,:)         = X(b,1:L,:) .* reshape(wL, [L 1 1]);
                X(b,end-L+1:end,:) = X(b,end-L+1:end,:) .* reshape(wR, [L 1 1]);
            end
        end
    end

    % === (3) Concatenate along time ===
    Y = reshape(permute(X, [2 1 3]), [B*T, nCh]);
    return
end

% -------------------------------------------------------------------------
% Mode 2: Default stitched alignment + polynomial detrend
% -------------------------------------------------------------------------
% Align the beginning of each block to the end of the preceding
% block, per channel. This preserves within-block shape and continuity
% across blocks better than forcing each block baseline to zero.
for b = 2:B
    prev_last  = squeeze(X(b-1,end,:));   % nCh x 1
    curr_first = squeeze(X(b,1,:));       % nCh x 1
    delta      = prev_last - curr_first;  % nCh x 1

    X(b,:,:) = X(b,:,:) + reshape(delta, [1 1 nCh]);
end

% Concatenate first, then remove the slow global drift that can accumulate
% from block-to-block offset matching.
Y = reshape(permute(X, [2 1 3]), [B*T, nCh]);
Y = local_poly_detrend(Y, opts.poly_order);

end


% =========================================================================
% Local helper
% =========================================================================
function Yd = local_poly_detrend(Y, poly_order)
% Remove a low-order polynomial trend per channel while preserving
% the fitted value at the first valid sample.

[Nt, nCh] = size(Y);
Yd = Y;
t  = (1:Nt)';

for ch = 1:nCh
    y  = Y(:,ch);
    ok = isfinite(y);

    if nnz(ok) <= poly_order + 1
        continue
    end

    p     = polyfit(t(ok), y(ok), poly_order);
    trend = polyval(p, t);

    first_idx = find(ok, 1, 'first');
    anchor    = trend(first_idx);

    Yd(:,ch) = y - trend + anchor;
end
end