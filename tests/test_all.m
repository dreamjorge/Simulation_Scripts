#!/usr/bin/env octave
% Test Suite Runner - Simulation_Scripts (CI & Local Wrapper)
% Uses portable_runner.m to ensure Octave 11 and MATLAB compatibility.

scriptPath = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptPath, '..', 'ParaxialBeams'));
addpath(scriptPath);

fprintf('=== Simulation_Scripts Test Suite (Portable) ===\n\n');

try
    status = portable_runner();
    
    if status == 0
        fprintf('\n=== ÉXITO: Todos los tests pasaron ===\n');
    else
        fprintf('\n=== FALLO: Se detectaron errores en la suite ===\n');
        % Exit with error if in Octave to signal CI failure
        if exist('OCTAVE_VERSION', 'builtin')
            exit(1);
        else
            error('Tests failed.');
        end
    end
catch ME
    fprintf('Error crítico ejecutando la suite: %s\n', ME.message);
    if exist('OCTAVE_VERSION', 'builtin')
        exit(1);
    else
        rethrow(ME);
    end
end

fprintf('\n=== Fin de la ejecución ===\n');
