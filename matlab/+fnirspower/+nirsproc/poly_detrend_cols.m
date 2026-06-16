function X = poly_detrend_cols(X, degree)
%NIRSPROC.POLY_DETREND_COLS  Column-wise polynomial detrend.
%  X = POLY_DETREND_COLS(X, degree)

t = (1:size(X,1))';
for k=1:size(X,2)
  pc = polyfit(t, X(:,k), degree);
  X(:,k) = X(:,k) - polyval(pc, t);
end
end
