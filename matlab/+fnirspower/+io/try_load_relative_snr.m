function val = try_load_relative_snr(glmdir)
val = NaN;
try
  sfile = dir(fullfile(glmdir,'*_Visual_Checkerboard_GLM_Relative_Beta_SNR.mat'));
  if ~isempty(sfile)
    [~,ix] = sort([sfile.datenum],'descend');
    S = load(fullfile(sfile(ix(1)).folder, sfile(ix(1)).name), 'relative_beta_SNR');
    if isfield(S,'relative_beta_SNR'), val = S.relative_beta_SNR; end
  end
catch
  val = NaN;
end
end