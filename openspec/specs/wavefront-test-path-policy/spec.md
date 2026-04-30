# Wavefront Test Path Policy Specification

## Purpose

Ensure Wavefront tests initialize paths from the actual repository root and do not emit false missing-path warnings.

## Requirements

### Requirement: Wavefront Test Uses Actual Repo Root

`tests/modern/test_Wavefront.m` MUST derive the repository root from its `tests/modern/` location as the directory two levels above the test file.

#### Scenario: Direct Wavefront invocation

- GIVEN `tests/modern/test_Wavefront.m` is run directly
- WHEN it initializes paths
- THEN it MUST add existing repo-root-based paths
- AND it MUST NOT attempt to add `tests/src/*` or `tests/ParaxialBeams/*` paths.

### Requirement: Wavefront Test Preserves Package Semantics

Wavefront test setup MUST preserve MATLAB/Octave package semantics for `+paraxial/` by adding the repository root rather than internal `+package` folders.

#### Scenario: Package parent path setup

- GIVEN `+paraxial/` exists at repository root
- WHEN Wavefront tests need canonical package classes
- THEN the repository root MUST be on the path
- AND internal `+paraxial` folders MUST NOT be added directly.

### Requirement: False Path Warnings Are Guarded

Repository guardrails SHOULD detect stale modern-test path setup that treats `tests/` as the repo root.

#### Scenario: Modern test path anti-pattern

- GIVEN a modern test computes `repoRoot = fullfile(testDir, '..')`
- WHEN guardrails scan modern test path setup
- THEN the guardrail SHOULD fail for tests that add project paths from that stale root
- AND it SHOULD pass when tests use `tests/modern/../..` or `setpaths()` semantics.
