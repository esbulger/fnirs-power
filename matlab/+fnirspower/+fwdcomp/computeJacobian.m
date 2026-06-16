function [J, DATA, nirsOut] = computeJacobian(nirsIn, modulationMHz)
%COMPUTEJACOBIAN Wrapper around NIRFASTer jacobian_FD.
% modulationMHz : scalar modulation frequency in MHz (0 for CW)

if nargin < 2, modulationMHz = 0; end

[J, DATA] = jacobian_FD(nirsIn, modulationMHz);
nirsOut = nirsIn; % return in case caller wants the mutated mesh
end