function init()
    % init - Initialize +paraxial package with convenient imports
    %
    % Usage:
    %   paraxial.init()  % Call once at start of session
    %
    % This adds all +paraxial/ subdirectories to the path and optionally
    % imports commonly used classes.
    %
    % Note: In MATLAB/Octave, adding +paraxial/ to path automatically
    % makes all subpackages accessible via import paraxial.* or
    % import paraxial.beams.* etc.

    scriptPath = fileparts(mfilename('fullpath'));

    % Add main package directory
    addpath(scriptPath);

    % Add all subpackages
    addpath(fullfile(scriptPath, '+beams'));
    addpath(fullfile(scriptPath, '+parameters'));
    addpath(fullfile(scriptPath, '+computation'));
    addpath(fullfile(scriptPath, '+propagation'));
    addpath(fullfile(scriptPath, '+propagation', '+field'));
    addpath(fullfile(scriptPath, '+propagation', '+rays'));
    addpath(fullfile(scriptPath, '+visualization'));

    fprintf('[+paraxial] Package initialized\n');
end