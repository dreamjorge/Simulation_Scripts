# paraxial-beam-factory Specification

## Purpose

BeamFactory creates beam instances from the canonical `+paraxial/` package namespace. The factory MUST support both `+paraxial/` canonical classes and legacy `src/` classes during the Strangler Fig migration transition, routing to `+paraxial/` by default while keeping `src/` as functional adapters with deprecation warnings.

## Requirements

### Requirement: Factory creates from +paraxial/ by default

BeamFactory.create() SHALL instantiate beam classes from `+paraxial/+beams/` when the class exists there. The factory SHALL use `which()` to resolve the canonical class location at runtime.

#### Scenario: Gaussian beam creation

- GIVEN BeamFactory and `+paraxial/+beams/GaussianBeam.m` exists and is loadable
- WHEN user calls `BeamFactory.create('gaussian', 100e-6, 632.8e-9)`
- THEN factory SHALL return an instance of `+paraxial/+beams/GaussianBeam`
- AND the instance SHALL be of class `GaussianBeam` in the `+paraxial` package

#### Scenario: Hermite beam creation

- GIVEN BeamFactory and `+paraxial/+beams/HermiteBeam.m` exists and is loadable
- WHEN user calls `BeamFactory.create('hermite', 100e-6, 632.8e-9, 'n', 2, 'm', 1)`
- THEN factory SHALL return an instance of `+paraxial/+beams/HermiteBeam`

### Requirement: Fallback to src/ during transition

If `+paraxial/` class is not on path, BeamFactory SHOULD fall back to `src/beams/` class. This supports incremental migration where only some beam classes have been migrated.

#### Scenario: Fallback when +paraxial/ class unavailable

- GIVEN `+paraxial/+beams/GaussianBeam.m` is NOT on the MATLAB/Octave path
- AND `src/beams/GaussianBeam.m` is on the path
- WHEN `BeamFactory.create('gaussian', 100e-6, 632.8e-9)` is called
- THEN factory SHALL fall back to `src/beams/GaussianBeam`
- AND emit a warning: "GaussianBeam from src/beams — migrate to +paraxial/"

### Requirement: Deprecation warnings on src/ classes

The `src/beams/` classes MUST emit a deprecation warning when instantiated directly (not via BeamFactory), guiding users to use BeamFactory or migrate to `+paraxial/`.

#### Scenario: Direct instantiation of src/ GaussianBeam

- GIVEN `src/beams/GaussianBeam.m` is on the path
- WHEN user directly instantiates `GaussianBeam(100e-6, 632.8e-9)` without using BeamFactory
- THEN the constructor SHALL emit warning: "src/beams/GaussianBeam is deprecated — use BeamFactory.create('gaussian', ...) or +paraxial/+beams/GaussianBeam directly"
- AND the instance SHALL still be returned (backward compatible)

### Requirement: Migration sequence verification

Each beam class migrated MUST have its tests pass in both Octave 11.1.0+ and MATLAB R2020b+ before the next class migration begins.

#### Scenario: GaussianBeam migration verification

- GIVEN `+paraxial/+beams/GaussianBeam.m` has been updated to canonical version
- WHEN `tests/test_all.m` is run
- THEN ALL GaussianBeam tests SHALL pass
- AND `tests/legacy_compat/run_legacy_compat.m` SHALL pass

### Requirement: No cross-namespace dependencies

`+paraxial/` beam classes MUST NOT import or depend on `src/` classes. All dependencies MUST be within `+paraxial/` or in shared utilities (`ParaxialBeams/`, `src/computation/`).

#### Scenario: Standalone +paraxial/ class verification

- GIVEN only `+paraxial/` package is on the path (not `src/beams/`)
- WHEN `BeamFactory.create('gaussian', 100e-6, 632.8e-9)` is called
- THEN it SHALL succeed
- AND the returned instance SHALL be fully functional

### Requirement: Factory beam type registry

BeamFactory SHALL maintain a mapping of beam type names to their canonical class names/locations, supporting the default public names returned by `BeamFactory.supportedTypes()`: 'gaussian', 'hermite', 'laguerre', 'elegant_hermite', 'elegant_laguerre', 'hankel', and 'hankel_hermite'. Public docs SHALL list only supported names as default API.

#### Scenario: Supported beam types

- GIVEN BeamFactory is properly initialized
- WHEN `BeamFactory.create()` is called with any supported type
- THEN factory SHALL resolve to correct `+paraxial/` class
- AND SHALL NOT throw "unrecognized beam type" error

#### Scenario: Documentation lists factory names

- GIVEN public docs list factory beam names
- WHEN the list is compared with the factory registry
- THEN every documented default name MUST be supported by BeamFactory
- AND unsupported historical aliases MUST be labeled legacy if mentioned

### Requirement: Public Documentation Matches Factory Canonical Role

Public docs MUST present `BeamFactory.create()` as the stable high-level constructor and `+paraxial/` as the canonical namespace for direct class access.

#### Scenario: Quick start uses stable entrypoint

- GIVEN a new user reads `README.md`
- WHEN they follow the primary beam creation example
- THEN the example SHOULD use `BeamFactory.create(...)`
- AND direct `src/beams` usage MUST NOT be shown as the recommended path

#### Scenario: Direct namespace usage is documented

- GIVEN documentation shows direct class construction
- WHEN the class belongs to a migrated beam
- THEN the example MAY use `paraxial.beams.<ClassName>(...)`
- AND SHOULD note Octave/MATLAB namespace limitations where relevant
