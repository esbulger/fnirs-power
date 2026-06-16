function [lambda_opt, lambda_vals, err_fit_vals, err_pen_vals, x_grid] = lcurve_lambda( ...
    obj_builder, lambda_vals, x0, lb, ub, fmin_opts)
%LCURVE_LAMBDA  Pick λ via L-curve corner using perpendicular distance.
%
%  [lambda_opt, lambda_vals, err_fit, err_pen, x_grid] = lcurve_lambda(obj_builder, lambda_vals, x0, lb, ub, fmin_opts)
%
%  obj_builder : @(lambda) returns [objFun, evalFun]
%                objFun(x) must return scalar objective used by fmincon
%                evalFun(x) must return [err, err_O, err_R, err_rat] (scalars)
%
%  err_fit is defined as err_O + wR*err_R to match the weighted objective.

  if nargin < 2 || isempty(lambda_vals), lambda_vals = logspace(-3,1,30); end
  if nargin < 3 || isempty(x0), x0 = [5.70; -1.91]; end
  if nargin < 4 || isempty(lb), lb = [0; -Inf]; end
  if nargin < 5 || isempty(ub), ub = [Inf; 0]; end
  if nargin < 6 || isempty(fmin_opts)
    fmin_opts = optimoptions('fmincon','Display','none','OptimalityTolerance',1e-8);
  end

  nL = numel(lambda_vals);
  err_fit_vals = zeros(nL,1);
  err_pen_vals = zeros(nL,1);
  x_grid = NaN(numel(x0), nL);

  for k = 1:nL
    lambda = lambda_vals(k);

    [objFun, evalFun] = obj_builder(lambda);

    x_opt = fmincon(objFun, x0, [],[],[],[], lb, ub, [], fmin_opts);
    x_grid(:,k) = x_opt;

    [~, err_O, err_R, err_rat] = evalFun(x_opt);

    % err_O/err_R/err_rat should be scalars; if not, fail loudly
    if ~isscalar(err_O) || ~isscalar(err_R) || ~isscalar(err_rat)
      error('evalFun must return scalar err_O, err_R, err_rat.');
    end

    % Weighted
    err_fit_vals(k) = err_O + err_R;
    err_pen_vals(k) = err_rat;
  end

  % Corner detection in log-log
  X = log10(max(eps, err_fit_vals));
  Y = log10(max(eps, err_pen_vals));
  p1 = [X(1) Y(1)];
  pN = [X(end) Y(end)];
  v  = pN - p1;
  v  = v / max(eps, norm(v));

  dist = zeros(nL,1);
  for i=1:nL
    w = [X(i) Y(i)] - p1;
    proj = p1 + dot(w,v)*v;
    dist(i) = norm([X(i) Y(i)] - proj);
  end

  [~, idx] = max(dist);
  lambda_opt = lambda_vals(idx);
end