function setpaths()
    % setpaths - Initialize path for Simulation_Scripts
    % Call this function before using the library, or add the
    % individual directories to your path.
    %
    % Modern structure (src/):
    %   addpath('src/beams');
    %   addpath('src/parameters');
    %   addpath('src/propagation/field');
    %   addpath('src/propagation/rays');
    %   addpath('src/visualization');
    %
    % Utilities (ParaxialBeams/):
    %   addpath('ParaxialBeams');
    %   addpath('ParaxialBeams/Addons');

    scriptPath = fileparts(mfilename('fullpath'));

    % Modern library structure (src/)
    addpath(fullfile(scriptPath, 'src', 'beams'));
    addpath(fullfile(scriptPath, 'src', 'parameters'));
    addpath(fullfile(scriptPath, 'src', 'propagation', 'field'));
    addpath(fullfile(scriptPath, 'src', 'propagation', 'rays'));
    addpath(fullfile(scriptPath, 'src', 'visualization'));

    % Utilities (retained in ParaxialBeams/)
    addpath(fullfile(scriptPath, 'ParaxialBeams'));
    addpath(fullfile(scriptPath, 'ParaxialBeams', 'Addons'));

    % Tests
    addpath(fullfile(scriptPath, 'tests'));

    fprintf('Path configured for Simulation_Scripts\n');
end
