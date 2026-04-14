% Compatible with GNU Octave and MATLAB
% Test Suite Runner - Simulation_Scripts (wrapper)
% Delegates to portable_runner() as the canonical implementation.

scriptPath = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptPath, '..', 'ParaxialBeams'));
addpath(scriptPath);

fprintf('=== Simulation_Scripts Test Suite (Portable Wrapper) ===\n\n');

try
    total_failed = portable_runner();

    if total_failed == 0
        fprintf('\n=== ÉXITO: Todos los tests pasaron ===\n');
    else
        fprintf('\n=== FALLO: Se detectaron %d errores en la suite ===\n', total_failed);
    end
catch ME
    fprintf('Error crítico ejecutando la suite: %s\n', ME.message);
    if ~exist('OCTAVE_VERSION', 'builtin')
        rethrow(ME);
    end
end

fprintf('\n=== Fin de la ejecución ===\n');
