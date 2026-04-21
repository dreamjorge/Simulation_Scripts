# Parameter Computation Layer Specification

## Purpose

Define the behavior of a centralized, stateless computation layer for beam parameter formulas so that parameter classes stop embedding duplicated numerical logic.

## Requirements

### Requirement: Provide stateless Gaussian core formulas

The system MUST provide static, stateless methods for Gaussian core quantities.

- GIVEN `w0`, `lambda`
- WHEN `rayleighDistance(w0, lambda)` is called
- THEN it MUST return `pi*w0^2/lambda`

- GIVEN `lambda`
- WHEN `waveNumber(lambda)` is called
- THEN it MUST return `2*pi/lambda`

- GIVEN `w0`, `z`, `zr`
- WHEN `waist(w0, z, lambda, zr)` is called
- THEN it MUST return `w0*sqrt(1+(z/zr)^2)`

- GIVEN `z=0`
- WHEN `radiusOfCurvature(0, zr)` is called
- THEN it MUST return `Inf`

#### Scenario: Waist at origin

- GIVEN `w0=100e-6`, `lambda=632.8e-9`, `z=0`
- WHEN waist is evaluated
- THEN result MUST equal `w0` within tolerance

#### Scenario: Waist at Rayleigh distance

- GIVEN `z=zr`
- WHEN waist is evaluated
- THEN result MUST equal `w0*sqrt(2)` within tolerance

### Requirement: Provide complex beam parameter helper

The system MUST provide complex beam parameter evaluation.

- GIVEN `z`, `zr`, `k`
- WHEN `complexBeamParameter(z, zr, k)` is called
- THEN it MUST return `1i*k/(2*(z+1i*zr))`

#### Scenario: Elegant helper equivalence

- GIVEN values from existing elegant parameter classes
- WHEN alpha is computed through delegation
- THEN alpha MUST match previous formula behavior within tolerance

### Requirement: Be framework-neutral and side-effect free

The computation layer MUST have no mutable object state and no side effects.

- GIVEN repeated calls with same inputs
- WHEN methods are evaluated
- THEN outputs MUST be deterministic and identical
