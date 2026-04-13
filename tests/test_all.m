#!/usr/bin/env octave
% Test Suite Runner - Simulation_Scripts
% Runs all modular tests

scriptPath = fileparts(mfilename('fullpath'));

fprintf('=== Simulation_Scripts Test Suite ===\n\n');

% Run modular tests
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

total_passed = 0;
total_failed = 0;

for i = 1:numel(testFiles)
    testFile = testFiles{i};
    fprintf('--- Running %s ---\n', testFile);
    try
        run(fullfile(scriptPath, testFile));
    catch ME
        fprintf('  FAIL: %s\n', ME.message);
    end
end

fprintf('\n=== All Tests Complete ===\n');
fprintf('Review output above for pass/fail counts per test file\n');
