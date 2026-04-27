# Pre-Merge Documentation Specification

## Purpose

This spec defines the documentation artifacts required to successfully merge `integration/pre-master` into `master`. The documentation MUST accurately reflect the current codebase structure and API contracts.

## Requirements

### Requirement: README.md Accuracy

The `README.md` file at the repository root SHALL reflect the actual repository structure. All file paths, class names, and examples MUST match existing files.

- GIVEN a developer reads `README.md`
- WHEN they follow the usage instructions
- THEN all commands MUST execute without "file not found" errors

#### Scenario: Usage Section Accuracy

- GIVEN `ParaxialBeams/` contains 27 `.m` files
- WHEN README describes `addpath ParaxialBeams`
- THEN this path MUST exist and be valid

#### Scenario: Class References Accuracy

- GIVEN `README.md` references `PhysicalConstants`, `GridUtils`, `FFTUtils`
- WHEN developer runs `help PhysicalConstants`
- THEN MATLAB/Octave MUST display class documentation

### Requirement: Architecture Documentation

The repository SHALL contain `docs/ARCHITECTURE.md` that documents:

- The class hierarchy (ParaxialBeam as abstract base)
- Strategy Pattern implementation (IPropagator interface)
- Factory Pattern (BeamFactory)
- Data flow: Grid → Beam → Propagator → Result

#### Scenario: Architecture Diagram

- GIVEN `docs/ARCHITECTURE.md` exists
- WHEN reader reviews the document
- THEN they can identify all 6 beam types and their relationships

### Requirement: Example Classification

Examples in `examples/` SHALL be classified as:

- `canonical`: Recommended for new users, fully tested
- `legacy`: Historical scripts, may contain deprecated patterns

#### Scenario: Canonical Example Identification

- GIVEN `examples/MainGauss_refactored.m` is marked canonical
- WHEN a new user runs the script
- THEN it MUST execute without errors AND produce expected beam output

## Acceptance Criteria

- [ ] `README.md` lists correct file paths matching `ls ParaxialBeams/*.m`
- [ ] `README.md` includes canonical example section
- [ ] `docs/ARCHITECTURE.md` exists with class diagrams
- [ ] `examples/` files have header comments indicating canonical/legacy status
