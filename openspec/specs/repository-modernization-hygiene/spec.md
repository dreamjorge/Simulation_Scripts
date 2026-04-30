# Repository Modernization Hygiene Specification

## Purpose

Define repository cleanup and documentation consistency rules for post-v2 modernization.

## Requirements

### Requirement: Single CI Source of Truth

The repository MUST document GitHub Actions as the canonical CI system when `.github/workflows/` contains active Octave, MATLAB, and release workflows.

#### Scenario: Stale CircleCI config exists

- GIVEN `.github/workflows/octave.yml`, `matlab.yml`, and `release.yml` exist
- WHEN repository hygiene is reviewed
- THEN `.circleci/config.yml` MUST be removed or explicitly documented as inactive
- AND public docs MUST NOT advertise CircleCI as canonical CI

### Requirement: Local Tooling Policy

Local agent/tooling metadata MUST be either intentionally versioned or ignored.

#### Scenario: `.opencode/` is present

- GIVEN `.opencode/package.json` exists in the working tree
- WHEN cleanup is applied
- THEN `.opencode/` MUST be added to `.gitignore` or committed with rationale
- AND the chosen policy MUST be documented in the cleanup change

### Requirement: Documentation Consistency

Public documentation MUST agree on supported runtime versions, canonical test commands, canonical API surface, and release packaging.

#### Scenario: Version support drift

- GIVEN `README.md` and `.atl/skill-registry.md` both mention Octave support
- WHEN docs are updated
- THEN both files MUST state the same supported Octave baseline
- AND any lower compatibility claim MUST be backed by tests or removed

#### Scenario: Test command drift

- GIVEN `tests/README.md` describes the canonical runner
- WHEN docs are updated
- THEN the command MUST match `tests/test_all.m` or `portable_runner()` usage in CI

### Requirement: Roadmap Reflects Current State

The project SHOULD keep a current roadmap separate from archived implementation plans.

#### Scenario: Historical plan remains at root

- GIVEN `plan.md` describes already-completed pre-merge work
- WHEN cleanup is applied
- THEN it SHOULD be archived or replaced by `docs/ROADMAP.md`
- AND active next steps SHOULD live in `docs/ROADMAP.md`
