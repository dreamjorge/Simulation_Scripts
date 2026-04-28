# matlab-toolbox-structure Specification

## Purpose

Define the structure for distributing Simulation_Scripts as a MATLAB toolbox (`.mltbx` bundle) installable via double-click or `matlab.addons.install`.

## Requirements

### Requirement: Toolbox Metadata

The system SHALL provide toolbox metadata in `toolbox_manifest.xml` (or equivalent) with fields: `name`, `version`, `author`, `email`, `description`, `website`.

#### Scenario: .mltbx contains valid metadata

- GIVEN a built `.mltbx` bundle
- WHEN `matlab.addons.install('simulation_scripts.mltbx')` is called
- THEN MATLAB SHALL read the toolbox name and version from the bundle
- AND display them in the Add-Ons Manager

### Requirement: Toolbox Installs All Namespaces

The system SHALL install `+paraxial/`, `ParaxialBeams/`, and `src/` folders into MATLAB's add-ons directory.

#### Scenario: Toolbox adds all namespaces to MATLAB path

- GIVEN `matlab.addons.install('simulation_scripts.mltbx')`
- WHEN installation completes
- THEN `+paraxial/+beams/GaussianBeam.m` SHALL be on MATLAB's path
- AND `ParaxialBeams/BeamFactory.m` SHALL be on MATLAB's path
- AND `src/beams/` SHALL be on MATLAB's path

#### Scenario: No test files in toolbox

- GIVEN the `.mltbx` is built
- THEN `tests/` folder SHALL NOT be included in the bundle
- AND `examples/` SHOULD be included for reference

### Requirement: Installation Script

The system SHALL provide `install.m` at toolbox root that MATLAB executes post-install.

#### Scenario: install.m runs after toolbox install

- GIVEN the toolbox is installed via double-click
- WHEN MATLAB starts after installation
- THEN `simulation_scripts_version()` SHALL be callable
- AND no manual `addpath` SHALL be required for basic usage

### Requirement: Uninstall Removes Paths

The system SHALL provide `uninstall.m` that MATLAB executes on toolbox uninstall.

#### Scenario: Uninstall removes all paths

- GIVEN `matlab.addons.uninstall('Simulation_Scripts')`
- WHEN uninstallation completes
- THEN `+paraxial` namespace SHALL NOT be on MATLAB's path
- AND all `src/` files SHALL be removed from path

### Requirement: Version Consistency

The system SHALL expose the same version string as the Git tag in both the toolbox metadata and `simulation_scripts_version()`.

#### Scenario: Version matches across package and code

- GIVEN Git tag `v2.0.0`
- WHEN the `.mltbx` is built and installed
- THEN `simulation_scripts_version()` SHALL return `'2.0.0'`
- AND the toolbox metadata version SHALL also be `'2.0.0'`