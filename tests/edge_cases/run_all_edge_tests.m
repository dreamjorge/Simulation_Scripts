% Compatible with GNU Octave and MATLAB
% Edge Case Test Suite Runner
%
% Runs all edge case tests for beam physics:
%   - GaussianBeam: z=0, origin, extreme parameters
%   - HankelHermite: z=0, r=0, theta=0, high orders
%   - HankelLaguerre: r=0, theta=0, z=0, l=0, p=0
%   - RayTracing: zero slope, negative slope, origin, extreme radial

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testDir = fileparts(mfilename('fullpath'));

fprintf('==============================================\n');
fprintf('   Edge Case Test Suite - Portable Runner\n');
fprintf('==============================================\n\n');

% Add paths
addpath(fullfile(repoRoot, 'ParaxialBeams'));
addpath(fullfile(repoRoot, 'src', 'beams'));
addpath(fullfile(repoRoot, 'src', 'parameters'));
addpath(fullfile(repoRoot, 'src', 'propagation', 'field'));
addpath(fullfile(repoRoot, 'src', 'propagation', 'rays'));
addpath(fullfile(repoRoot, 'src', 'visualization'));
addpath(fullfile(repoRoot, 'legacy', 'compat'));

% List of edge case test files
edgeTests = {
    'test_GaussianBeam_edge.m',
    'test_HankelHermite_edge.m',
    'test_HankelLaguerre_edge.m',
    'test_RayTracing_extreme.m'
};

totalPassed = 0;
totalFailed = 0;

for i = 1:numel(edgeTests)
    testFile = fullfile(testDir, edgeTests{i});
    [~, testName, ~] = fileparts(edgeTests{i});
    
    fprintf('----------------------------------------------\n');
    fprintf('Running: %s\n', testName);
    fprintf('----------------------------------------------\n');
    
    try
        run(testFile);
        
        % Parse output for passed/failed counts
        % This is a simple heuristic - tests print their own summary
        fprintf('\n');
        
    catch ME
        fprintf('  [ERROR] Test crashed: %s\n', ME.message);
        totalFailed = totalFailed + 1;
    end
end

fprintf('==============================================\n');
fprintf('   Edge Case Suite Summary\n');
fprintf('==============================================\n');
fprintf('Test files executed: %d\n', numel(edgeTests));
fprintf('Note: See individual test output above for pass/fail counts\n');
fprintf('==============================================\n');