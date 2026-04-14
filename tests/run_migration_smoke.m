% Compatible with GNU Octave and MATLAB
% Migration smoke runner for canonical examples + portable tests.

scriptPath = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptPath, '..');

examplesToRun = {
    fullfile(repoRoot, 'examples', 'canonical', 'MainGauss_refactored.m')
    fullfile(repoRoot, 'examples', 'canonical', 'MainMultiMode.m')
};

fprintf('=== Migration Smoke Runner ===\n\n');

exampleFailures = 0;
for i = 1:numel(examplesToRun)
    examplePath = examplesToRun{i};
    [~, exampleName] = fileparts(examplePath);
    fprintf('--- Running example: %s ---\n', exampleName);
    prevDir = pwd;
    try
        cd(repoRoot);
        run(examplePath);
        cd(prevDir);
        fprintf('  PASS: %s\n', exampleName);
    catch ME
        cd(prevDir);
        fprintf('  FAIL: %s (%s)\n', exampleName, ME.message);
        exampleFailures = exampleFailures + 1;
    end
    close all;
end

fprintf('\n--- Running portable test suite ---\n');
prevDir = pwd;
cd(scriptPath);
testFailures = portable_runner();
cd(prevDir);

totalFailures = exampleFailures + testFailures;

fprintf('\n=== Smoke Summary ===\n');
fprintf('Example failures: %d\n', exampleFailures);
fprintf('Test failures: %d\n', testFailures);

if totalFailures ~= 0
    error('Migration smoke failed (%d total failures).', totalFailures);
end

fprintf('Migration smoke passed.\n');
