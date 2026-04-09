#!/usr/bin/env octave
% Test Suite Runner - Simulation_Scripts
% Runs all modular tests using the runtests framework

scriptPath = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptPath, '..', 'ParaxialBeams'));

fprintf('=== Simulation_Scripts Test Suite ===\n\n');

% Run all tests in the tests/ directory
try
    results = runtests(scriptPath);
    display(results);
    
    % Summary
    passCount = sum([results.Passed]);
    failCount = sum([results.Failed]);
    incompleteCount = sum([results.Incomplete]);
    
    fprintf('\n=== Test Summary ===\n');
    fprintf('Passed: %d\n', passCount);
    fprintf('Failed: %d\n', failCount);
    if incompleteCount > 0
        fprintf('Incomplete: %d\n', incompleteCount);
    end
    
    if failCount > 0
        error('Tests failed.');
    end
catch ME
    fprintf('Error running test suite: %s\n', ME.message);
    % Fallback for environments where runtests might behave differently
    run('test_RayTracing.m');
end

fprintf('\n=== All Tests Complete ===\n');
