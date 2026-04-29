# Delta for CI Test Status

## MODIFIED Requirements

### Requirement: Portable Runner Failure Propagation

CI jobs that execute `portable_runner()` MUST fail when the returned failure count is non-zero. Documentation that references CI MUST identify GitHub Actions as the active CI system and MUST describe runner failure propagation consistently with workflow commands.
(Previously: CI failure propagation was specified, but documentation alignment with active workflow files was not required.)

#### Scenario: Portable runner succeeds

- GIVEN `portable_runner()` returns `0`
- WHEN the Octave portable CI step runs
- THEN the step MUST complete successfully
- AND `portable-tests.log` SHOULD still be written

#### Scenario: Portable runner reports failures

- GIVEN `portable_runner()` returns a value greater than `0`
- WHEN the Octave portable CI step runs
- THEN the step MUST raise an error
- AND the GitHub Actions job MUST fail
- AND `portable-tests.log` SHOULD still be uploaded by the artifact step

#### Scenario: Runner crashes before returning status

- GIVEN Octave raises an exception while executing the runner
- WHEN the workflow command is piped through `tee`
- THEN the pipeline MUST fail through shell `pipefail`

#### Scenario: Documentation references CI

- GIVEN public docs describe CI or test execution
- WHEN cleanup is applied
- THEN docs MUST match `.github/workflows/octave.yml` and `.github/workflows/matlab.yml`
- AND stale CI systems MUST NOT be described as active
