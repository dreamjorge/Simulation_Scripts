# version-reporting Specification

## Purpose

Provide a version-reporting function so installed packages and toolboxes can report their semantic version programmatically.

## Requirements

### Requirement: Version Function Exists

The system SHALL provide `simulation_scripts_version()` function accessible after installation.

#### Scenario: Function is callable after install

- GIVEN Simulation_Scripts package/toolbox is installed
- WHEN user calls `v = simulation_scripts_version()`
- THEN `v` SHALL be a character vector
- AND it SHALL NOT error

### Requirement: Version Format

The returned version SHALL be a valid semantic version string: `MAJOR.MINOR.PATCH` (e.g., `'2.0.0'`).

#### Scenario: Returns valid semver

- GIVEN tag `v2.0.0` is applied
- WHEN `simulation_scripts_version()` is called
- THEN the returned string SHALL match the regex `^\d+\.\d+\.\d+$`
- AND it SHALL be `'2.0.0'`

### Requirement: Version Source is Git Tag

The version SHALL be derived from the most recent Git tag on the current commit.

#### Scenario: Version matches git tag

- GIVEN `git describe --tags --abbrev=0` returns `v2.1.3`
- WHEN `simulation_scripts_version()` is called
- THEN it SHALL return `'2.1.3'`

#### Scenario: No tag returns fallback

- GIVEN the commit has no tags
- WHEN `simulation_scripts_version()` is called
- THEN it SHALL return `'0.0.0'` or `'unreleased'`

### Requirement: Implementation Location

The version function SHALL be implemented in `+paraxial/init.m`.

#### Scenario: Version lives in paraxial init

- GIVEN the package is installed
- WHEN `which('simulation_scripts_version')` is called
- THEN it SHALL point to `+paraxial/init.m`