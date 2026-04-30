function setpaths()
    % setpaths - Initialize path for Simulation_Scripts
    % Call this function before using the library, or add the
    % individual directories to your path.
    %
    % DUAL-PATH SUPPORT:
    % This function adds both the legacy 'src/' paths and the modern
    % '+paraxial/' package paths. Users can choose which to use.
    %
    % Legacy structure (src/):
    %   addpath('src/beams');
    %   addpath('src/parameters');
    %   ...
    %
    % Modern package (+paraxial/):
    %   addpath(repoRoot);  % package parent enables 'import paraxial.*'
    %   import paraxial.beams.GaussianBeam
    %
    % Utilities (ParaxialBeams/):
    %   addpath('ParaxialBeams');
    %   addpath('ParaxialBeams/Addons');

    scriptPath = fileparts(mfilename('fullpath'));

    %% Legacy library structure (src/) — backward compatibility
    addpath(fullfile(scriptPath, 'src', 'beams'));
    addpath(fullfile(scriptPath, 'src', 'parameters'));
    addpath(fullfile(scriptPath, 'src', 'computation'));
    addpath(fullfile(scriptPath, 'src', 'propagation', 'field'));
    addpath(fullfile(scriptPath, 'src', 'propagation', 'rays'));
    addpath(fullfile(scriptPath, 'src', 'visualization'));

    %% +paraxial namespace package
    % MATLAB/Octave packages are resolved by adding the parent directory.
    % Do not add internal +package folders directly; doing so pollutes
    % unqualified class resolution and produces Octave package-path warnings.
    addpath(scriptPath);
    addpath(fullfile(scriptPath, 'ParaxialBeams'));     % BeamFactory and utilities
    addpath(fullfile(scriptPath, 'ParaxialBeams', 'Addons'));

    %% Legacy compatibility aliases
    addpath(fullfile(scriptPath, 'legacy', 'compat'));

    %% Tests
    addpath(fullfile(scriptPath, 'tests'));

    fprintf('Path configured for Simulation_Scripts (dual-path: src/ + +paraxial/)\n');
end
