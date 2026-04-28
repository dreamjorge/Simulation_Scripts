# octave-package-structure Specification

## Purpose

Define the structure and metadata for distributing Simulation_Scripts as an Octave package via `.tar.gz` archive, enabling `pkg install` without manual path management.

## Requirements

### Requirement: Package Metadata

The system SHALL provide `DESCR.ini` at package root with fields: `Name`, `Version`, `Author`, `Maintainer`, `Description`, `Website`, `License`, `Depends`.

#### Scenario: DESCR.ini exists at package root

- GIVEN a built `.tar.gz` package
- WHEN the archive is inspected
- THEN `DESCR.ini` SHALL exist at the root of the archive
- AND it SHALL contain all required fields

#### Scenario: Version matches git tag

- GIVEN `DESCR.ini` declares `Version = 2.0.0`
- WHEN the package is built from commit tagged `v2.0.0`
- THEN the version SHALL match exactly

### Requirement: Namespace Preservation

The system SHALL preserve the full `+paraxial/` namespace structure within the package archive.

#### Scenario: Package installs +paraxial namespace correctly

- GIVEN `pkg install simulation_scripts-2.0.0.tar.gz`
- WHEN Octave starts
- THEN `+paraxial/+beams/GaussianBeam.m` SHALL be loadable via `which`
- AND `BeamFactory.create` SHALL route to `+paraxial/+beams/`

#### Scenario: ParaxialBeams utilities on path

- GIVEN the package is installed
- THEN `BeamFactory`, `GridUtils`, `PhysicalConstants` SHALL be on Octave's load path
- AND `tests/` folder SHALL NOT be included in the package

### Requirement: Installation Script

The system SHOULD provide an `install.m` script at package root that runs automatically after `pkg install`.

#### Scenario: install.m executes on package install

- GIVEN `pkg install simulation_scripts-2.0.0.tar.gz`
- WHEN the package installation completes
- THEN Octave SHALL automatically execute `install.m` if present
- AND paths SHALL be added without manual `addpath` calls

### Requirement: Uninstall Script

The system SHALL provide an `uninstall.m` script that removes all added paths.

#### Scenario: Uninstall removes paths cleanly

- GIVEN `pkg uninstall simulation_scripts`
- WHEN uninstallation completes
- THEN `+paraxial` namespace SHALL NOT be on Octave's path
- AND no artifacts SHALL remain in Octave's package directory

### Requirement: Package Build Reproducibility

The system SHALL produce identical `.tar.gz` output for the same Git commit regardless of build machine.

#### Scenario: Build is deterministic

- GIVEN commit `abc123`
- WHEN package is built on machine A and machine B
- THEN both resulting `.tar.gz` files SHALL be byte-identical
- OR if non-deterministic (timestamps), a `CHECKSUMS.md5` file SHALL be included