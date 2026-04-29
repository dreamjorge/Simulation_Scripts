# Delta for paraxial-beam-factory

## ADDED Requirements

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

## MODIFIED Requirements

### Requirement: Factory beam type registry

BeamFactory SHALL maintain a mapping of beam type names to their canonical class names/locations, supporting: 'gaussian', 'hermite', 'laguerre', 'elegant-hermite', 'elegant-laguerre', 'hankel-hermite', 'hankel-laguerre'. Public docs SHALL list only supported names as default API.
(Previously: the registry requirement did not require public docs to stay aligned with supported names.)

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
