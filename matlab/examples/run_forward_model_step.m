%% RUN_FORWARD_MODEL_STEP  (combined + cleaned)
%  - Loads head mesh + NIRx montage (Standard_probeInfo.mat)
%  - Runs NIRFASTer forward model + saves outputs
%  - Uses fnirspower.paths() for consistent project paths
import fnirspower.*

P = fnirspower.setup_paths();

%% ---------------------- User inputs ----------------------
montageName   = "temporal_nq";          % folder name under montages/ (or layouts/)
lambdas_nm    = [760 850];        % wavelengths
modulationMHz = 0;                % CW
ri            = 1.4;

% Head model assets
meshMatPath     = fullfile(P.icbm_mesh_dir,'ICBM_mesh_5layer.mat');        % contains genmesh.node/genmesh.ele
baseNirfastMesh = fullfile(P.icbm_mesh_dir,'nirs_mesh_ICBM_5layer');       % NIRFASTer mesh prefix (no ext)

% NIRx probe info
% (choose ONE of these patterns depending on where you store them)
probeInfoMat = fullfile(P.montages, char(montageName), 'Standard_probeInfo.mat');
if exist(probeInfoMat,'file')~=2
  error('Could not find Standard_probeInfo.mat for montage "%s".', montageName);
end

% For testing with built-in NIRx montage
% probeInfoMat = "C:\Users\elibu\Documents\NIRx\Configurations\Montages\Headband_8x8\Standard_probeInfo.mat";

% Extinction coefficients file for NIRFASTer
excoefPath = fullfile(P.thirdparty,'NIRFASTer-master','toolbox','common','excoef.txt');
if exist(excoefPath,'file')~=2
  error('Missing excoef.txt at: %s', excoefPath);
end

% Output locations
outPrefix = fullfile(P.forward, sprintf('nirs_mesh_ICBM_5layer_%s', montageName));  % model bundle prefix
save_dir  = fullfile(P.derivatives,'forward');

% Optode placement (only if your montage coordinates need scaling/offset onto the head)
optPlacement = struct( ...
  'scale',       10, ...
  'offset',      [0 18 0], ...
  'zero_center', false);

debugPlots = true;
%% ---------------------- Run ----------------------
out = fnirspower.pipeline.run_forward_model( ...
    'meshMatPath', meshMatPath, ...
    'baseNirfastMesh', baseNirfastMesh, ...
    'probeInfoMat', probeInfoMat, ...
    'outPrefix', outPrefix, ...
    'excoefPath', excoefPath, ...
    'lambdas', lambdas_nm, ...
    'modulationMHz', modulationMHz, ...
    'ri', ri, ...
    'debugPlots', debugPlots, ...
    'optPlacement', optPlacement);


%% Save model
fprintf("Saving model to %s...", out.outMeshPath)

nirs = out.nirs;
J_760 = out.J_760;
J_850 = out.J_850;
DATA760 = out.DATA760;
DATA850 = out.DATA850;
epsilon = out.epsilon;
Jspec = out.Jspec;

save([out.outMeshPath '_nirsmodel.mat'], ...
    'nirs', 'J_760', 'J_850', 'DATA760', 'DATA850', 'epsilon', 'Jspec', '-v7.3');