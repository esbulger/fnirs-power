function [raw, Fs, time, Y, link, wl_idx, t1_idx, t2_idx] = load_trim_snirf(snirf_path, buffers)
%LOAD_TRIM_SNIRF  Load SNIRF, compute trim window from stimuli, and slice.
%  [raw, Fs, time, Y, link, wl_idx, t1_idx, t2_idx] = ...
%      nirsproc.load_trim_snirf(snirf_path, buffers)
%
%  Inputs
%    snirf_path : path to a single .snirf file
%    buffers    : [pre post] seconds to include around first/last onset
%                 default = [20 30]
%
%  Outputs
%    raw     : struct returned by loadsnirf (with .nirs.* fields)
%    Fs      : sampling rate (Hz)
%    time    : trimmed time vector (starts at 0)
%    Y       : trimmed dataTimeSeries (T×M)
%    link    : M/2 × 2 matrix [src det] for channel set (wavelength-agnostic)
%    wl_idx  : M × 1 vector, wavelength index (1 for 760, 2 for 850, ...)
%    t1_idx  : start index within the original full-length time series
%    t2_idx  : end index within the original full-length time series
%
%  Notes
%   - Stimuli can contain multiple conditions; this function uses the
%     earliest first onset and latest last onset across all conditions to
%     define the trimming window, then re-bases each condition to t=0.
%   - If a condition has no events inside the window after trimming, it is
%     kept but with times shifted (users can ignore downstream).
%
%  Requirements: loadsnirf (JSNIRFy/EasyH5 compatible)

arguments
  snirf_path (1,:) char
  buffers (1,2) double = [20 30]
end

% --- Load SNIRF ---
raw = loadsnirf(snirf_path);

% Ensure time is a column vector
full_time = raw.nirs.data.time(:);
Fs = 1 / (full_time(2) - full_time(1));

% --- Determine trim window from stimuli ---
if ~isfield(raw.nirs, 'stim') || isempty(raw.nirs.stim)
  error('load_trim_snirf:NoStim', 'No stimulus info found in %s', snirf_path);
end

% Collect first/last onsets across all conditions
first_on = inf; last_on = -inf;
for k = 1:numel(raw.nirs.stim)
  t_on = raw.nirs.stim(k).data(:,1);
  if ~isempty(t_on)
    first_on = min(first_on, t_on(1));
    last_on  = max(last_on,  t_on(end));
  end
end
if ~isfinite(first_on) || ~isfinite(last_on)
  error('load_trim_snirf:EmptyStim', 'Stim fields present but empty onsets.');
end

pre  = buffers(1);
post = buffers(2);

% Boundaries in seconds
t1_sec = max(0, first_on - pre);
t2_sec = min(last_on + post, full_time(end));

% Convert to indices using time vector (robust to non-zero start)
t1_idx = find(full_time >= t1_sec, 1, 'first');
t2_idx = find(full_time <= t2_sec, 1, 'last');

if isempty(t1_idx) || isempty(t2_idx) || t2_idx <= t1_idx
  error('load_trim_snirf:BadWindow', 'Invalid trim window [%g %g] sec.', t1_sec, t2_sec);
end

% --- Slice data/time and rebase to 0 ---
Y = raw.nirs.data.dataTimeSeries(t1_idx:t2_idx, :);
time = full_time(t1_idx:t2_idx);
time = time - time(1);  % start at 0 s

% --- Shift stimulus times to match the new origin ---
dt0 = full_time(t1_idx);  % original seconds at new t=0
for k = 1:numel(raw.nirs.stim)
  if ~isempty(raw.nirs.stim(k).data)
    raw.nirs.stim(k).data(:,1) = raw.nirs.stim(k).data(:,1) - dt0;
  end
end

% --- Build link matrix and wavelength indices ---
mlist = raw.nirs.data.measurementList;
% Make sure we can handle struct array or cell array of structs
if iscell(mlist); mlist = [mlist{:}]; end
src = [mlist.sourceIndex]';
det = [mlist.detectorIndex]';
wl_idx = [mlist.wavelengthIndex]';

% Return an M/2×2 link (per channel, wavelength-agnostic) by taking the
% first wavelength occurrences; we keep wl_idx separately for split.
% Assume same ordering per-ch: 1..nCh for wl1, then 1..nCh for wl2
nMeas = numel(wl_idx);
% Heuristic: nCh = min(count(wl==1), count(wl==2), ...)
wl_vals = unique(wl_idx);
counts = arrayfun(@(w) sum(wl_idx==w), wl_vals);
nCh = min(counts);
link = [src(1:nCh) det(1:nCh)];

end