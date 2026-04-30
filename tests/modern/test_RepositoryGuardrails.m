% Guardrail: keep repository cleanup and canonical API documentation aligned

repoRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(repoRoot, 'ParaxialBeams'));

fprintf('=== Repository Guardrail Tests ===\n\n');
passed = 0;
failed = 0;

% -------------------------------------------------------------------------
% Local tooling and CI hygiene
% -------------------------------------------------------------------------
gitignorePath = fullfile(repoRoot, '.gitignore');
circleCiPath = fullfile(repoRoot, '.circleci', 'config.yml');
if exist(gitignorePath, 'file')
    gitignoreContent = fileread(gitignorePath);
else
    gitignoreContent = '';
end

if ~isempty(strfind(gitignoreContent, '.opencode/'))
    fprintf('  PASS: .opencode/ is ignored as local tooling\n');
    passed = passed + 1;
else
    fprintf('  FAIL: .opencode/ is not listed in .gitignore\n');
    failed = failed + 1;
end

if ~exist(circleCiPath, 'file')
    fprintf('  PASS: stale CircleCI config is absent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: stale CircleCI config still exists at %s\n', circleCiPath);
    failed = failed + 1;
end

% -------------------------------------------------------------------------
% Documentation invariants
% -------------------------------------------------------------------------
readmePath = fullfile(repoRoot, 'README.md');
registryPath = fullfile(repoRoot, '.atl', 'skill-registry.md');
if exist(readmePath, 'file')
    readmeContent = fileread(readmePath);
else
    readmeContent = '';
end
if exist(registryPath, 'file')
    registryContent = fileread(registryPath);
else
    registryContent = '';
end

if ~isempty(strfind(readmeContent, 'GitHub Actions is the canonical CI system'))
    fprintf('  PASS: README documents GitHub Actions as canonical CI\n');
    passed = passed + 1;
else
    fprintf('  FAIL: README does not document GitHub Actions as canonical CI\n');
    failed = failed + 1;
end

if ~isempty(strfind(readmeContent, '+paraxial/')) && ~isempty(strfind(readmeContent, 'deprecated'))
    fprintf('  PASS: README documents canonical +paraxial/ and deprecated src/ surfaces\n');
    passed = passed + 1;
else
    fprintf('  FAIL: README canonical/deprecated namespace documentation is incomplete\n');
    failed = failed + 1;
end

if ~isempty(strfind(registryContent, 'Octave 11.1.0+')) && isempty(strfind(registryContent, 'Octave 6+'))
    fprintf('  PASS: skill registry Octave baseline matches public docs\n');
    passed = passed + 1;
else
    fprintf('  FAIL: skill registry Octave baseline is stale or inconsistent\n');
    failed = failed + 1;
end

% -------------------------------------------------------------------------
% Canonical examples must not regress to direct src/beams usage
% -------------------------------------------------------------------------
canonicalDir = fullfile(repoRoot, 'examples', 'canonical');
exampleFiles = dir(fullfile(canonicalDir, '*.m'));
exampleViolations = {};

for i = 1:numel(exampleFiles)
    filePath = fullfile(canonicalDir, exampleFiles(i).name);
    content = fileread(filePath);

    if ~isempty(strfind(content, "addpath('src/beams')")) || ...
       ~isempty(strfind(content, 'addpath("src/beams")'))
        exampleViolations{end+1} = filePath; %#ok<AGROW>
    end
end

if isempty(exampleViolations)
    fprintf('  PASS: canonical examples avoid direct src/beams path usage\n');
    passed = passed + 1;
else
    fprintf('  FAIL: canonical examples use direct src/beams paths\n');
    for i = 1:numel(exampleViolations)
        fprintf('    - %s\n', exampleViolations{i});
    end
    failed = failed + 1;
end

% -------------------------------------------------------------------------
% BeamFactory supported registry stays explicit
% -------------------------------------------------------------------------
expectedTypes = {'gaussian', 'hermite', 'laguerre', ...
                 'elegant_hermite', 'elegant_laguerre', ...
                 'hankel', 'hankel_hermite'};
supportedTypes = BeamFactory.supportedTypes();
missingTypes = {};

for i = 1:numel(expectedTypes)
    if ~any(strcmp(supportedTypes, expectedTypes{i}))
        missingTypes{end+1} = expectedTypes{i}; %#ok<AGROW>
    end
end

if isempty(missingTypes)
    fprintf('  PASS: BeamFactory supported type registry contains expected public names\n');
    passed = passed + 1;
else
    fprintf('  FAIL: BeamFactory supported type registry is missing expected names\n');
    for i = 1:numel(missingTypes)
        fprintf('    - %s\n', missingTypes{i});
    end
    failed = failed + 1;
end

fprintf('\n=== Repository Guardrails: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Repository guardrails failed: %d violation group(s).', failed);
end
