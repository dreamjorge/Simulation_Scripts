# Beam API Contract Specification

## Purpose

This spec documents the formal contract that ALL ParaxialBeam subclasses MUST implement. This contract ensures propagators (FFT, Ray Trace, Analytic) work with any beam type without knowing internal details.

## Requirements

### Requirement: Unified opticalField Interface

ALL beam types MUST implement `opticalField(obj, X, Y, z)` method that returns the complex field on a Cartesian grid.

- GIVEN any ParaxialBeam subclass (GaussianBeam, HermiteBeam, LaguerreBeam, etc.)
- WHEN `opticalField(X, Y, z)` is called with 2D matrices X, Y and scalar z
- THEN it MUST return a [Ny x Nx] complex array representing the optical field

#### Scenario: GaussianBeam opticalField

- GIVEN `beam = GaussianBeam(100e-6, 632.8e-9)`
- WHEN `field = beam.opticalField(X, Y, 0)`
- THEN `field` MUST be [Ny x Nx] complex array where Gaussian amplitude peaks at origin

#### Scenario: HermiteBeam opticalField

- GIVEN `beam = HermiteBeam(100e-6, 632.8e-9, 1, 0)`
- WHEN `field = beam.opticalField(X, Y, 0)`
- THEN `field` MUST exhibit Hermite-Gaussian mode shape

### Requirement: Dynamic Parameter Evaluation

ALL beam types MUST implement `getParameters(obj, z)` that returns a GaussianParameters object evaluated at position z.

- GIVEN any ParaxialBeam subclass
- WHEN `getParameters(z)` is called with a given z
- THEN it MUST return an object with at least: Waist, RayleighDistance, GouyPhase, RadiusOfCurvature

#### Scenario: Gaussian Parameters at z=0

- GIVEN `beam = GaussianBeam(100e-6, 632.8e-9)`
- WHEN `params = beam.getParameters(0)`
- THEN `params.Waist` MUST equal `100e-6`

#### Scenario: Gaussian Parameters at z≠0

- GIVEN `beam = GaussianBeam(100e-6, 632.8e-9)`
- WHEN `params = beam.getParameters(0.05)`
- THEN `params.Waist` MUST be > `100e-6` (beam diverges)

### Requirement: Beam Name Identification

ALL beam types MUST implement `beamName(obj)` that returns a string identifier.

- GIVEN any ParaxialBeam subclass
- WHEN `beamName()` is called
- THEN it MUST return a string like 'gaussian', 'hermite_3_2', 'laguerre_1_0'

#### Scenario: HermiteBeam Name with Mode Indices

- GIVEN `beam = HermiteBeam(w0, lambda, 3, 2)`
- WHEN `name = beam.beamName()`
- THEN `name` MUST be 'hermite_3_2'

### Requirement: Shared State

ALL beam types MUST store `Lambda` (wavelength in meters) and `k` (wave number = 2π/λ).

- GIVEN any ParaxialBeam subclass
- WHEN instantiated with lambda parameter
- THEN `obj.Lambda == lambda` AND `obj.k == 2*pi/lambda`

## Acceptance Criteria

- [ ] `GaussianBeam` implements all 3 interface methods
- [ ] `HermiteBeam` implements all 3 interface methods
- [ ] `LaguerreBeam` implements all 3 interface methods
- [ ] `ElegantHermiteBeam` implements all 3 interface methods
- [ ] `ElegantLaguerreBeam` implements all 3 interface methods
- [ ] `HankelLaguerre` implements all 3 interface methods
- [ ] All propagators accept any ParaxialBeam subclass via IPropagator interface
