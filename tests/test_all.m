#!/usr/bin/env octave
% Test Suite Runner - Simulation_Scripts (CI & Local Wrapper)
% Uses portable_runner.m to ensure Octave 11 and MATLAB compatibility.

scriptPath = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptPath, '..', 'ParaxialBeams'));
addpath(scriptPath);

fprintf('=== Simulation_Scripts Test Suite (Portable) ===\n\n');


% modular test list from utility-classes
testFiles = {
    'test_PhysicalConstants.m'
    'test_GridUtils.m'
    'test_FFTUtils.m'
    'test_GaussianParameters.m'
    'test_HermiteParameters.m'
    'test_LaguerreParameters.m'
    'test_ElegantHermiteParameters.m'
    'test_ElegantLaguerreParameters.m'
    'test_GaussianBeam.m'
    'test_HermiteBeam.m'
    'test_LaguerreBeam.m'
    'test_ElegantHermiteBeam.m'
    'test_ElegantLaguerreBeam.m'
    'test_HankelLaguerre.m'
    'test_CylindricalRay.m'
    'test_OpticalRay.m'
    'test_AnalysisUtils.m'
};

try
    % Uses the modular runner logic
    scriptPath = fileparts(mfilename('fullpath'));
    addpath(scriptPath);
    
    total_failed = 0;
    for i = 1:numel(testFiles)
        testFile = testFiles{i};
        fprintf('--- Running %s ---\n', testFile);
        try
            run(fullfile(scriptPath, testFile));
        catch ME
            fprintf('  FAIL: %s\n', ME.message);
            total_failed = total_failed + 1;
        end
    end
    
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
