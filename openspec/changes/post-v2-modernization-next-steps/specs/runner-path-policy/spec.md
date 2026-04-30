# Runner Path Policy Specification

## Purpose

Define how development and test runners expose canonical, deprecated, utility, and legacy paths during the Strangler Fig migration.

## Requirements

### Requirement: Canonical Package Parent Path

The test runner MUST add the repository root, not internal `+package` folders, so `+paraxial/` resolves through MATLAB/Octave package semantics.

#### Scenario: Portable runner initializes canonical package

- GIVEN `+paraxial/` exists at repo root
- WHEN `tests/portable_runner.m` initializes paths
- THEN repo root MUST be added to the path
- AND internal `+paraxial` subfolders MUST NOT be added directly.

### Requirement: Deprecated Paths Are Explicit

Deprecated `src/` paths MAY remain in the runner for compatibility, but they MUST be labeled as deprecated/transitional compatibility paths.

#### Scenario: Runner keeps compatibility paths

- GIVEN legacy compatibility tests still require `src/` adapters
- WHEN the runner adds `src/beams` or related folders
- THEN comments MUST identify them as deprecated compatibility paths
- AND comments MUST NOT call them the modern library paths.

### Requirement: setpaths Remains Preferred Dev Setup

Developer-facing docs SHOULD prefer `setpaths()` for local setup and direct test invocation.

#### Scenario: Docs describe manual setup

- GIVEN docs show local MATLAB/Octave path setup
- WHEN users follow the recommended setup
- THEN `setpaths()` SHOULD be the first recommendation
- AND package parent path semantics SHOULD be preserved.
