function P = paths()
%FNIRSPOWER.PATHS Return centralized, project-relative paths.
%
%  P = fnirspower.paths()
%
%  Description
%  -----------
%  This function defines the recommended repository and workspace layout
%  used by the example scripts and convenience workflows. All paths are
%  resolved relative to the location of this file:
%
%      <repo>/matlab/+fnirspower/paths.m
%
%  The returned struct includes paths to:
%    - repository and MATLAB roots
%    - workspace assets (models, layouts, montages, rawdata)
%    - output directories (derivatives)
%    - bundled third-party toolboxes
%    - fnirspower package submodules
%
%  Notes
%  -----
%  - This function is intended mainly for example scripts and default
%    project-layout resolution.
%  - Core pipeline functions should still accept explicit file paths and
%    should not require this layout.
%
%  Returns
%  -------
%  P : struct of project-relative paths
%
%  Example
%  -------
%  P = fnirspower.paths();
%  mesh_file = fullfile(P.icbm_mesh_dir, 'ICBM_mesh_5layer.mat');
%
%  See also
%  --------
%  fnirspower.pipeline.run_forward_model

here      = mfilename('fullpath');     % this file
pkg_dir   = fileparts(here);           % .../+fnirspower
matlabdir = fileparts(pkg_dir);        % .../matlab
repo_root = fileparts(matlabdir);      % .../<repo root>

% ---- Project roots
P.root       = repo_root;
P.matlab     = matlabdir;

% ---- Workspace assets
P.workspace     = fullfile(repo_root,'workspace');
P.models        = fullfile(P.workspace,'models');
P.icbm_mesh_dir = fullfile(P.models,'icbm2009cna','mesh');      % << meshes
P.icbm_segment  = fullfile(P.models,'icbm2009cna','segment');
P.forward       = fullfile(P.models,'icbm2009cna','forward');
P.montages      = fullfile(P.workspace,'montages');

% ---- Data
P.rawdata    = fullfile(P.workspace,'rawdata');                  % raw inputs
P.examples   = fullfile(repo_root,'examples');

% ---- Layouts / figures / docs
P.layouts    = fullfile(P.workspace,'layouts');
P.docs       = fullfile(repo_root,'docs');

% ---- Outputs
P.derivatives = fullfile(repo_root,'workspace','derivatives');   % recommended
P.figures    = fullfile(P.derivatives,'figures');

% ---- MATLAB / package folders
P.thirdparty = fullfile(matlabdir,'thirdparty');
P.easyh5     = fullfile(P.thirdparty,'EasyH5');
P.jsnirfy    = fullfile(P.thirdparty,'jsnirfy-master');
P.fromNIRS   = fullfile(P.thirdparty,'fromNIRSToolbox');
P.iso2mesh   = fullfile(P.thirdparty,'iso2mesh');
P.fieldtrip  = fullfile(P.thirdparty,'Fieldtrip');
P.nirfast  = fullfile(P.thirdparty,'NIRFASTer-master');

% ---- MATLAB package folders
P.fnirspower = fullfile(matlabdir,'+fnirspower');
P.nirsproc   = fullfile(P.fnirspower,'+nirsproc');
P.glmproc    = fullfile(P.fnirspower,'+glmproc');
P.recon      = fullfile(P.fnirspower,'+recon');
P.variance   = fullfile(P.fnirspower,'+variance');
P.absmag     = fullfile(P.fnirspower,'+absmag');
P.fwdcomp    = fullfile(P.fnirspower,'+fwdcomp');
P.pipeline   = fullfile(P.fnirspower,'+pipeline');
P.helpers    = fullfile(P.fnirspower,'+helpers');
P.io         = fullfile(P.fnirspower,'+io');
P.measpred   = fullfile(P.fnirspower,'+measpred');

% ---- Resolve any ".." and normalize
fn = fieldnames(P);
for k = 1:numel(fn)
    P.(fn{k}) = char(java.io.File(P.(fn{k})).getCanonicalPath());
end
end