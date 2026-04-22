% Guardrail: prevent deprecated Hankele* aliases in canonical/modern surfaces

repoRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));

fprintf('=== Legacy Alias Guardrail Tests ===\n\n');
passed = 0;
failed = 0;

forbiddenAliases = {'HankeleHermite', 'HankeleLaguerre'};

targets = {
    fullfile(repoRoot, 'examples', 'canonical'),
    fullfile(repoRoot, 'tests', 'modern')
};

violations = {};

for t = 1:numel(targets)
    dirPath = targets{t};
    files = dir(fullfile(dirPath, '*.m'));

    for i = 1:numel(files)
        filePath = fullfile(dirPath, files(i).name);

        % Ignore this guardrail test file to avoid self-matching constants.
        if strcmp(files(i).name, 'test_LegacyAliasGuardrail.m')
            continue;
        end

        content = fileread(filePath);

        for a = 1:numel(forbiddenAliases)
            alias = forbiddenAliases{a};
            % Match token usage to reduce false positives.
            hasAlias = ~isempty(regexp(content, ['\<' alias '\>'], 'once'));

            if hasAlias
                violations{end+1} = sprintf('%s -> %s', filePath, alias); %#ok<AGROW>
            end
        end
    end
end

if isempty(violations)
    fprintf('  PASS: no deprecated Hankele* aliases in canonical/modern surfaces\n');
    passed = passed + 1;
else
    fprintf('  FAIL: found deprecated Hankele* alias usage\n');
    for i = 1:numel(violations)
        fprintf('    - %s\n', violations{i});
    end
    failed = failed + 1;
end

fprintf('\n=== Legacy Alias Guardrail: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Legacy alias guardrail failed: %d violation group(s).', failed);
end
