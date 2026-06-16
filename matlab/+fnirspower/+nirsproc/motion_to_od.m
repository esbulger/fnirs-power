function dOD = motion_to_od(Y, Fs)
%NIRSPROC.MOTION_TO_OD  Apply TDDR motion correction and convert to ΔOD.
%  dOD = MOTION_TO_OD(Y, Fs)

Ycorr = TDDR(Y, Fs);
dOD = -bsxfun(@minus, real(log(Ycorr)), log(mean(Ycorr)));
end
