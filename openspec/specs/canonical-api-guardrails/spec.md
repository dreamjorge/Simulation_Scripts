# Canonical API Guardrails Specification

## Purpose

Define structural checks that keep new project surfaces aligned with `+paraxial/` and `BeamFactory` as canonical APIs.

## Requirements

### Requirement: Canonical Examples Avoid Deprecated Direct Imports

Canonical examples MUST NOT directly add or instantiate deprecated `src/` beam paths.

#### Scenario: Canonical example path audit

- GIVEN files under `examples/canonical/*.m`
- WHEN guardrail tests scan their contents
- THEN examples MUST NOT call `addpath('src/beams')`
- AND SHOULD use `setpaths()`, `BeamFactory.create(...)`, or `paraxial.*` APIs

### Requirement: Documentation Names Canonical API Surface

Public documentation MUST describe `+paraxial/` and `BeamFactory` as canonical for new code.

#### Scenario: README API audit

- GIVEN `README.md` is reviewed by a guardrail test or manual audit
- WHEN canonical API sections are inspected
- THEN `+paraxial/` MUST be named as canonical namespace
- AND `src/` MUST be marked deprecated or transitional

### Requirement: Factory Remains Canonical Entrypoint

BeamFactory documentation and tests MUST keep supported beam type names explicit.

#### Scenario: Beam type registry audit

- GIVEN `ParaxialBeams/BeamFactory.m` supports beam creation
- WHEN guardrails inspect expected public beam names
- THEN the supported set MUST include gaussian, hermite, laguerre, elegant_hermite, elegant_laguerre, hankel, and hankel_hermite
- AND docs MUST NOT recommend unsupported aliases as default API

### Requirement: Deprecated Surfaces Stay Isolated

Deprecated compatibility surfaces MAY remain, but new docs and examples MUST NOT present them as default usage.

#### Scenario: Legacy example boundary

- GIVEN files under `examples/legacy/`
- WHEN public docs mention them
- THEN they MUST be described as archive, generator, research, or compatibility material
- AND they MUST NOT be the first recommended path for new users
