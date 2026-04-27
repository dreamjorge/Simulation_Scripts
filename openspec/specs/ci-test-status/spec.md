# CI Test Status Specification

## Purpose

Define how CI jobs propagate MATLAB/Octave test runner status.

## Requirements

### Requirement: Portable Runner Failure Propagation

CI jobs that execute `portable_runner()` MUST fail when the returned failure count is non-zero.

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
