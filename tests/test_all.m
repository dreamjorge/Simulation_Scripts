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
