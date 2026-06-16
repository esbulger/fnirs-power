function [HR, bad, nCh] = estimate_hr_snr_bad(dOD_hp, Fs, snr_thr)
%NIRSPROC.ESTIMATE_HR_SNR_BAD  HR (Hz) from spectrum and SNR-based bad channels.
%  [HR, bad, nCh] = ESTIMATE_HR_SNR_BAD(dOD_hp, Fs, snr_thr)

% channel count per wavelength
nCh = size(dOD_hp,2)/2;

[freq,amp] = get_spec(Fs, mean(dOD_hp(:, std(dOD_hp)<0.08), 2));
amp(freq<0.8 | freq>1.5) = 0; [~,iHR] = max(amp); HR = freq(iHR);

snr = zeros(nCh,1);
for jj=1:nCh
  [f2,a2] = get_spec(Fs, dOD_hp(:, nCh+jj));
  sPow = sum(a2(f2>HR-0.1 & f2<HR+0.1));
  nPow = sum(a2(f2>1.5*HR-0.1 & f2<1.5*HR+0.1));
  snr(jj) = sPow/max(eps,nPow);
end
bad = find(snr<snr_thr);
end
