# Parameter Model Delegation Specification

## Purpose

Ensure existing parameter classes preserve public API behavior while delegating numerical calculations to the new stateless computation layer.

## Requirements

### Requirement: Preserve parameter class public API

All existing constructors, public methods, and dependent properties in parameter classes MUST remain available.

- GIVEN existing consumer code
- WHEN calling `GaussianParameters`, `HermiteParameters`, `LaguerreParameters`, `ElegantHermiteParameters`, `ElegantLaguerreParameters`
- THEN no signature changes or removals are allowed

#### Scenario: Gaussian snapshot API compatibility

- GIVEN `params = GaussianParameters(z, w0, lambda)`
- WHEN reading `params.Waist`, `params.GouyPhase`, `params.Radius`
- THEN values MUST be available and numerically equivalent to pre-refactor behavior

#### Scenario: Gaussian dynamic API compatibility

- GIVEN `params = GaussianParameters(z, w0, lambda)`
- WHEN calling `params.waist(z2)`, `params.gouyPhase(z2)`, `params.radius(z2)`
- THEN results MUST remain numerically equivalent to pre-refactor behavior

### Requirement: Delegate computations from parameter classes

Parameter classes MUST not own duplicated numerical formulas once the refactor is complete; they MUST delegate formula evaluation to the computation layer.

- GIVEN parameter family methods (Gaussian/Hermite/Laguerre/Elegant)
- WHEN computing waist, phase, curvature, or alpha quantities
- THEN the values MUST be obtained via computation-layer methods

#### Scenario: Hermite modal phase delegation

- GIVEN `HermiteParameters(n,m)`
- WHEN `phiPhase(z)` is called
- THEN formula source MUST be delegated and output remain `(n+m)*psi(z)`

#### Scenario: Laguerre modal phase delegation

- GIVEN `LaguerreParameters(l,p)`
- WHEN `phiPhase(z)` is called
- THEN output MUST remain `(|l|+2p)*psi(z)`

### Requirement: Preserve legacy Hermite helper entrypoint

The historical `HermiteParameters.getHermiteSolutions(...)` entrypoint MUST remain callable.

- GIVEN legacy script calls to `HermiteParameters.getHermiteSolutions(nu, x)`
- WHEN evaluated after refactor
- THEN call MUST succeed and return outputs equivalent to previous behavior

#### Scenario: Shim delegation behavior

- GIVEN both old and new helper entrypoints
- WHEN called with same inputs
- THEN outputs MUST match within tolerance
