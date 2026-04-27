# Legacy Policy Specification

## Purpose

Define how legacy examples are classified, maintained, and eventually removed without disrupting current research workflows.

## Requirements

### Requirement: Example Classification

The system SHALL classify all examples into exactly one of three categories:

- **canonical**: Recommended entrypoints for new users. Use `BeamFactory.create()` and `setpaths()`.
- **legacy/archive**: Old API reference, retained for historical reference. May be removed in future major versions.
- **legacy/research**: Thesis or paper-specific scripts with hardcoded parameters. Preserved for reproducibility.
- **legacy/generators**: Scripts that generate specific figures or data for papers. Depend on external data inputs.

### Requirement: Canonical Examples Contract

Canonical examples MUST:

- Use `BeamFactory.create()` for beam instantiation
- Use `setpaths()` for path initialization
- Be runnable in Octave 11.1.0+ and MATLAB R2020b+
- Have verified output against theory

### Requirement: Legacy Retention Policy

Legacy examples in `legacy/archive/`:

- SHALL NOT be actively maintained
- SHALL be kept for historical reference
- MAY be removed in future major versions (v3.0.0+)

Legacy examples in `legacy/research/` and `legacy/generators/`:

- SHALL be preserved for reproducibility of published results
- SHALL NOT be actively maintained
- MAY be archived to external storage in future major versions

### Requirement: Legacy Migration Path

When upgrading legacy scripts to canonical:

1. Replace `addpath('ParaxialBeams')` with `setpaths()`
2. Replace old constructor patterns with `BeamFactory.create()`
3. Replace legacy `addpath` chains with `setpaths()`

## Scenarios

### Scenario: Classify a new example

- GIVEN a new example script
- WHEN the maintainer evaluates it for inclusion
- THEN they MUST classify it as canonical, legacy/archive, legacy/research, or legacy/generators

### Scenario: Deprecate a legacy example

- GIVEN an example in `legacy/archive/`
- WHEN it is no longer needed for reference
- THEN it MAY be removed in the next major version release with a deprecation notice

### Scenario: Migrate a legacy script to canonical

- GIVEN `examples/legacy/archive/MainGauss.m` using old API
- WHEN a maintainer updates it to canonical
- THEN the updated script MUST use `BeamFactory.create()` and `setpaths()` and be moved to `examples/canonical/`
