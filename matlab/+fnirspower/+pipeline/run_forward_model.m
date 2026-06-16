function out = run_forward_model(opts)
%FNIRSPOWER.PIPELINE.RUN_FORWARD_MODEL
% Build a two-wavelength NIRFAST/NIRFASTer forward model.
%
% This function loads a tetrahedral head mesh and NIRx probe geometry,
% projects the optodes onto the head surface, constructs a NIRFAST mesh,
% computes wavelength-specific Jacobians at 760 and 850 nm, and converts
% them to a spectral HbO/HbR Jacobian.
%
% Usage
% -----
% out = fnirspower.pipeline.run_forward_model( ...
%     meshMatPath='/path/to/head_mesh.mat', ...
%     baseNirfastMesh='/path/to/base_nirfast_mesh', ...
%     probeInfoMat='/path/to/Standard_probeInfo.mat', ...
%     outPrefix='/path/to/output/nirs_mesh_montage', ...
%     excoefPath='/path/to/excoef.txt');
%
% Required name-value arguments
% -----------------------------
% meshMatPath :
%     Path to a MAT file containing the tetrahedral head mesh used for
%     optode projection. The file must contain the mesh variables expected
%     by fnirspower.fwdcomp.loadHeadTetraMesh.
%
% baseNirfastMesh :
%     Path or prefix identifying the base NIRFAST/NIRFASTer mesh used by
%     fnirspower.fwdcomp.buildNirsMesh.
%
% probeInfoMat :
%     Path to a MAT file containing a variable named probeInfo.
%
% outPrefix :
%     Output prefix used when saving the optode-placed NIRFAST mesh.
%
% excoefPath :
%     Path to the extinction-coefficient text file used to construct the
%     spectral HbO/HbR Jacobian.
%
% Optional name-value arguments
% -----------------------------
% lambdas :
%     Wavelengths in nanometres. The current implementation requires
%     exactly [760 850]. Default: [760 850].
%
% modulationMHz :
%     Modulation frequency passed to the NIRFAST Jacobian calculation, in
%     megahertz. Default: 0.
%
% ri :
%     Refractive index passed to buildNirsMesh. Default: 1.4.
%
% debugPlots :
%     If true, plot the source and detector coordinates from the reloaded
%     NIRFAST mesh. Default: false.
%
% optPlacement :
%     Scalar struct passed to
%     fnirspower.fwdcomp.placeAndProjectOptodes. It may contain placement
%     options such as coordinate scaling, offsets, and centering behavior.
%     Default: struct().
%
% fwhm :
%     Optional smoothing parameter passed to buildNirsMesh. Default: [].
%
% Output
% ------
% out :
%     Struct containing:
%
%       out.nirs
%           Reloaded NIRFAST/NIRFASTer mesh containing the placed optodes.
%
%       out.J_760
%           Jacobian calculated at 760 nm.
%
%       out.J_850
%           Jacobian calculated at 850 nm.
%
%       out.DATA760
%           Forward-model output calculated at 760 nm.
%
%       out.DATA850
%           Forward-model output calculated at 850 nm.
%
%       out.epsilon
%           Extinction-coefficient matrix for the requested wavelengths.
%
%       out.Jspec
%           Spectral Jacobian expressed in HbO/HbR space.
%
%       out.outMeshPath
%           Prefix used to save the optode-placed NIRFAST mesh.
%
% Files written
% -------------
% buildNirsMesh writes the optode-placed NIRFAST mesh using outPrefix.
% This function returns the remaining forward-model results in out but does
% not independently save the out struct.
%
% Notes
% -----
% The optical properties are currently fixed to the package defaults for
% 760 and 850 nm and assume five tissue regions:
%
%   1 = white matter
%   2 = grey matter
%   3 = cerebrospinal fluid
%   4 = skull
%   5 = scalp
%
% FieldTrip is not required, but the NIRFAST/NIRFASTer functions and the
% fnirspower forward-model helpers must be available on the MATLAB path.
%
% See also:
%   fnirspower.fwdcomp.loadHeadTetraMesh
%   fnirspower.fwdcomp.placeAndProjectOptodes
%   fnirspower.fwdcomp.buildNirsMesh
%   fnirspower.fwdcomp.computeSpectralJacobian


arguments
    opts.meshMatPath (1,:) char
    opts.baseNirfastMesh (1,:) char
    opts.probeInfoMat (1,:) char
    opts.outPrefix (1,:) char
    opts.excoefPath (1,:) char
    opts.lambdas (1,:) double = [760 850]
    opts.modulationMHz (1,1) double = 0
    opts.ri (1,1) double = 1.4
    opts.debugPlots (1,1) logical = false
    opts.optPlacement struct = struct()
    opts.fwhm double = []
end

import fnirspower.fwdcomp.*

% 1) Load head mesh (tetra) for projection
head = loadHeadTetraMesh(opts.meshMatPath);

% 2) Load probeInfo
S = load(opts.probeInfoMat, 'probeInfo');
if ~isfield(S, 'probeInfo')
    error('probeInfoMat must contain a variable named "probeInfo".');
end
probeInfo = S.probeInfo;

% add optode placement parameters into higher level 'opts'
opts.scale = opts.optPlacement.scale;
opts.offset = opts.optPlacement.offset;
opts.zero_center = opts.optPlacement.zero_center;

% 3) Place & project optodes
[srcProj, detProj, link] = placeAndProjectOptodes(head, probeInfo, opts);

% 4) Build NIRFASTer mesh with optodes and save
disp("Defining geometry and building mesh... This may take several minutes...")
outMeshPath = opts.outPrefix;
nirsPlaced = buildNirsMesh( ...
    opts.baseNirfastMesh, srcProj, detProj, link, outMeshPath, ...
    struct('ri', opts.ri, 'fwhm', opts.fwhm));

% 5) Reload to ensure saved mesh coordinates are used
nirs = load_mesh(outMeshPath);

% 6) Optional quick plot of optode locations - NIRFAST will slighyly change
% locations
if opts.debugPlots
    hold on
    scatter3(nirs.source.coord(:,1), nirs.source.coord(:,2), nirs.source.coord(:,3), 50, 'filled');
    scatter3(nirs.meas.coord(:,1),   nirs.meas.coord(:,2),   nirs.meas.coord(:,3),   50, 'filled');
    axis equal
    view(90,0)
    legend('S', 'D')
    hold off
end

% 7) Optical properties per wavelength (Eggebrecht-like defaults)
% Regions: 1=WM, 2=GM, 3=CSF, 4=Skull, 5=Scalp

props760.r(1).mua = 0.0167; props760.r(1).mus = 1.1908;
props760.r(2).mua = 0.0180; props760.r(2).mus = 0.8359;
props760.r(3).mua = 0.0040; props760.r(3).mus = 0.30;
props760.r(4).mua = 0.0116; props760.r(4).mus = 0.94;
props760.r(5).mua = 0.0170; props760.r(5).mus = 0.74;

props850.r(1).mua = 0.0208; props850.r(1).mus = 1.0107;
props850.r(2).mua = 0.0192; props850.r(2).mus = 0.6726;
props850.r(3).mua = 0.0040; props850.r(3).mus = 0.30;
props850.r(4).mua = 0.0139; props850.r(4).mus = 0.84;
props850.r(5).mua = 0.0190; props850.r(5).mus = 0.64;

% 8) Compute Jacobians explicitly for 760 and 850 nm
disp("Computing Jacobians. This may take several minutes...")

if numel(opts.lambdas) ~= 2 || ~all(sort(opts.lambdas) == [760 850])
    error('run_forward_model currently expects opts.lambdas = [760 850].');
end

nirs760 = setOpticalProps(nirs, props760);
[J_760, DATA760] = computeJacobian(nirs760, opts.modulationMHz);

nirs850 = setOpticalProps(nirs, props850);
[J_850, DATA850] = computeJacobian(nirs850, opts.modulationMHz);

% 9) Spectral Jacobian (HbO/HbR)
epsilon = fnirspower.fwdcomp.readExtinctionCoeffs(opts.excoefPath, opts.lambdas);
Jspec = fnirspower.fwdcomp.computeSpectralJacobian({J_760; J_850}, epsilon);

% Return
out = struct();
out.nirs = nirs;
out.J_760 = J_760;
out.J_850 = J_850;
out.DATA760 = DATA760;
out.DATA850 = DATA850;
out.epsilon = epsilon;
out.Jspec = Jspec;
out.outMeshPath = outMeshPath;
end