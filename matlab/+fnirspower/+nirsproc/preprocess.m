function S = preprocess(snirf_path, model, params)
%NIRSPROC.PREPROCESS  Motion→OD→filters→badch→DPF→Hb + stim vector (modular).
%  S = NIRSPROC.PREPROCESS(snirf_path, model, params)
%
%  See also: nirsproc.load_trim_snirf, nirsproc.motion_to_od,
%            nirsproc.poly_detrend_cols, nirsproc.attenuate_spikes,
%            nirsproc.filter_signals, nirsproc.estimate_hr_snr_bad,
%            nirsproc.split_wavelengths, nirsproc.sd_distances,
%            nirsproc.compute_dpfs, nirsproc.od_to_hb, nirsproc.build_stimvec
%
%  Author: Eli Bulger (modular refactor)
%  Version: 0.2.0

arguments
  snirf_path (1,:) char
  model struct
  params.buffers (1,2) double = [20 30]
  params.hp_hz double = 0.01
  params.bp_hz (1,2) double = [0.01 0.2]
  params.snr_thr double = 1.5
  params.poly_detrend double = 2
end

import fnirspower.*

fprintf('Preprocessing\n')

fprintf('\t Loading data...\n')
% 1) Load & trim
[raw, Fs, time, Y, link, wl_idx, t1_idx, t2_idx] = nirsproc.load_trim_snirf(snirf_path, params.buffers);

fprintf('\t Applying motion correction (TDDR) and converting to OD...\n')
% 2) Motion correction → OD
[dOD] = nirsproc.motion_to_od(Y, Fs);

% 3) Polynomial detrend (optional)
if params.poly_detrend > 0
    fprintf('\t Applying polynomial detrending...\n')
    dOD = nirsproc.poly_detrend_cols(dOD, params.poly_detrend);
end

% 4) Spike attenuation
fprintf('\t Attentuating spikes...\n')
[dOD] = nirsproc.attenuate_spikes(dOD, 30, 30, 0.02);

% 5) Filtering (HP and BP variants)
fprintf('\t Filtering...\n')
[dOD_hp, dOD_bp] = nirsproc.filter_signals(dOD, Fs, params.hp_hz, params.bp_hz);

% 6) Estimate HR + SNR → bad channels
fprintf('\t Estimating bad channels...\n')
[HR, bad, nCh] = nirsproc.estimate_hr_snr_bad(dOD_hp, Fs, params.snr_thr);

% 7) Split wavelengths (760/850)
[dOD_760, dOD_850] = nirsproc.split_wavelengths(dOD_hp, wl_idx);

% 8) Source–detector distances
[dist_ch] = nirsproc.sd_distances(raw, link);

% 9) Differential pathlengths (DPF) from Jacobians
[DPF_760, DPF_850] = nirsproc.compute_dpfs( ...
    model, model.epsilon_mm_uM, model.brain_nodes_idx, wl_idx, dist_ch, nCh);

% 10) Convert ΔOD → HbO/HbR
fprintf('\t Converting to HbX via mBLL...\n')
C = nirsproc.coeffs_defaults();
[HBO, HBD] = nirsproc.od_to_hb(dOD_hp, nCh, DPF_760, DPF_850, dist_ch, C);

% 11) Stimulus vector
fprintf('\t Building stim vector...\n')
stimvec = nirsproc.build_stimvec(raw, Fs, numel(time));
% ------------------------------------------------------------
% 12) Accelerometer regressor (real if available, else proxy)
% ------------------------------------------------------------
fprintf('\t Building accelerometer regressor...\n');

[snirf_dir, snirf_base, ~] = fileparts(snirf_path);
nirs_path = fullfile(snirf_dir, [snirf_base '.nirs']);

has_real_acc = false; acc_reg = []; acc_meta = struct();

if exist(nirs_path,'file') == 2
  try
    Saux = load(nirs_path, '-mat', 'aux');   % aux dims: [T_full x 2 x (axes)] or similar
    if isfield(Saux,'aux') && ~isempty(Saux.aux)
      % Combine axes, band-pass
      [b_acc, a_acc] = butter(3, [0.01, 0.2]/(Fs/2));
      acc = [squeeze(Saux.aux(:,1,:)), squeeze(Saux.aux(:,2,:))];
      % Prefer exact trim by index if available
      if exist('t1_idx','var') && ~isempty(t1_idx) && ~isempty(t2_idx)
        acc_win = acc(t1_idx:t2_idx, :);
      else
        % Fallback: match length by simple crop/pad (rare path)
        Twant = numel(time);
        Tfull = size(acc,1);
        a = max(1, floor((Tfull - Twant)/2) + 1);
        z = min(Tfull, a + Twant - 1);
        acc_win = acc(a:z, :);
        if size(acc_win,1) ~= Twant
          % pad to match, last value repeat
          acc_win(end+1:Twant, :) = repmat(acc_win(end,:), [Twant - size(acc_win,1), 1]);
        end
      end
      acc_filt = filtfilt(b_acc, a_acc, acc_win);
      acc_reg  = mean(acc_filt, 2);       % average axes → single regressor
      has_real_acc = true;
      acc_meta.path = nirs_path;
      acc_meta.filter = [0.01 0.2];
      acc_meta.trim_idx = [exist('t1_idx','var') && ~isempty(t1_idx), exist('t2_idx','var') && ~isempty(t2_idx)];
    end
  catch ME
    warning('Failed to load/use aux from %s: %s', nirs_path, ME.message);
  end
end

if ~has_real_acc
  fprintf('\t  No accelerometer found → proxy from mean |dOD_hp|\n');
  acc_reg = mean(dOD_hp, 2);
  acc_meta.path = '';
  acc_meta.filter = [params.bp_hz(1) params.bp_hz(2)];
  acc_meta.proxy = true;
else
  fprintf('\t  Accelerometer found and used\n');
  acc_meta.proxy = false;
end

% ------------------------------------------------------------
% Package
% ------------------------------------------------------------
S = struct('Fs', Fs, 'time', time, 'dOD_hp', dOD_hp, 'dOD_bp', dOD_bp, ...
           'dOD_760', dOD_760, 'dOD_850', dOD_850, 'HBO', HBO, 'HBD', HBD, ...
           'bad_channels', bad, 'nCh', nCh, 'stimvec', stimvec, ...
           'dist_ch', dist_ch, 'DPF_760', DPF_760, 'DPF_850', DPF_850, ...
           'HR_Hz', HR, ...
           'acc_reg', acc_reg, 'has_real_acc', has_real_acc, 'acc_meta', acc_meta);

fprintf('Preprocessing complete\n');

end
