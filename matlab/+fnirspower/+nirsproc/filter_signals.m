function [dOD_hp, dOD_bp] = filter_signals(dOD, Fs, hp_hz, bp_hz)
%NIRSPROC.FILTER_SIGNALS  Produce high-pass and band-pass filtered versions.
%  [dOD_hp, dOD_bp] = FILTER_SIGNALS(dOD, Fs, hp_hz, bp_hz)

[b,a] = butter(3, hp_hz/(Fs/2), 'high');   dOD_hp = filtfilt(b, a, dOD);
[b,a] = butter(3, bp_hz/(Fs/2));            dOD_bp = filtfilt(b, a, dOD);
end
