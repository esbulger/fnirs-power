function P = setup_paths()
%FNIRSPOWER.SETUP_PATHS Add fnirspower and third-party dependencies to the MATLAB path.
%
%  P = fnirspower.setup_paths()
%
%  Adds:
%    - fnirspower MATLAB root
%    - selected third-party dependencies
%
%  Notes
%  -----
%  - FieldTrip is added only at the top-level folder. Call ft_defaults
%    separately if needed.
%  - This function is intended as a convenience for examples and interactive use.

P = fnirspower.paths();

% Main MATLAB and third party code roots
addpath(P.matlab);
addpath(P.thirdparty);

if exist(P.workspace, 'dir') == 7
    addpath(genpath(P.workspace));
end
% Third-party packages that should include subfolders
if exist(P.easyh5, 'dir') == 7
    addpath(genpath(P.easyh5));
end

if exist(P.jsnirfy, 'dir') == 7
    addpath(genpath(P.jsnirfy));
end

if exist(P.fromNIRS, 'dir') == 7
    addpath(genpath(P.fromNIRS));
end

% Third-party packages added at top level only
if exist(P.iso2mesh, 'dir') == 7
    addpath(genpath(P.iso2mesh));
end

% Third-party packages added at top level only
if exist(P.nirfast, 'dir') == 7
    addpath(genpath(P.nirfast));
end

if exist(P.fieldtrip, 'dir') == 7
    addpath(P.fieldtrip);
    ft_defaults;
end
end