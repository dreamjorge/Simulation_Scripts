# Delta for Runner Path Policy

## ADDED Requirements

### Requirement: Direct Modern Test Invocation Uses Repo Root

Modern test scripts that set up project paths for direct invocation SHOULD derive the actual repository root before adding project directories.

#### Scenario: Modern test direct setup

- GIVEN a test under `tests/modern/` is run without `portable_runner()`
- WHEN it adds project paths directly
- THEN it SHOULD derive repo root as two levels above the test file
- AND it SHOULD avoid adding non-existent `tests/src/*` or `tests/ParaxialBeams/*` paths.
