function x = invert_woodbury(A, y, bad, alpha)
%RECON.INVERT_WOODBURY  Tikhonov recon with spatially varying penalty.
%  x = RECON.INVERT_WOODBURY(A, y, bad, alpha)
%
%  Inputs
%    A     : channels × vertices Jacobian
%    y     : channels × 1 vector (e.g., task beta per channel)
%    bad   : vector of bad-channel indices to drop
%    alpha : scalar Tikhonov parameter
%
%  Output
%    x     : vertices × 1 reconstructed map
%
%  Notes
%    Uses LL = diag(sqrt(sum(A.^2) + 0.1)) as in original script and
%    Woodbury identity for efficiency.

arguments
  A double
  y (:,1) double
  bad double = []
  alpha (1,1) double = 0.01
end

A(bad,:) = []; y(bad,:) = [];

LL = spdiags( ...
    sqrt(sum(A.^2)' + 0.1), ...
    0, size(A,2), size(A,2));

LLi = inv(LL);
At = A * LLi;

x = LLi * ( ...
    At' / alpha * y ...
    - At' * ( ...
        (inv(alpha * eye(size(A,1)) + At * At') * At) ...
        * (At' / alpha * y) ...
    ) ...
);

end
