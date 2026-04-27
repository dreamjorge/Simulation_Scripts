# Delta: Deprecation Propagation

## MODIFIED Requirements

### Requirement: Deprecation Warning Emission

Classes in `src/` MUST emit a deprecation warning when instantiated, indicating the user should use `+paraxial/` or `BeamFactory.create()` instead.

(Previously: Only `src/beams/` classes emitted deprecation warnings)

#### Scenario: Deprecated propagation class instantiation

- GIVEN `src/propagation/field/FFTPropagator(grid, lambda)` is called
- WHEN the class is instantiated
- THEN a deprecation warning SHALL be emitted: "src/propagation/field/FFTPropagator is deprecated. Use BeamFactory.create() or +paraxial/+propagation/+field/FFTPropagator directly."

#### Scenario: Deprecated ray propagation class instantiation

- GIVEN `src/propagation/rays/RayTracePropagator(grid)` is called
- WHEN the class is instantiated
- THEN a deprecation warning SHALL be emitted with the same pattern

### Requirement: Deprecation Warning Content

Deprecation warnings MUST include:

1. The deprecated path (e.g., `src/propagation/field/FFTPropagator`)
2. The recommended replacement path (e.g., `+paraxial/+propagation/+field/FFTPropagator` or `BeamFactory.create()`)
3. A warning ID for programmatic filtering (e.g., `BeamFactory:deprecatedSrc`)

## ADDED Requirements

### Requirement: Propagation Class Deprecation Coverage

The following classes MUST emit deprecation warnings:

- `src/propagation/field/IPropagator.m`
- `src/propagation/field/FFTPropagator.m`
- `src/propagation/field/AnalyticPropagator.m`
- `src/propagation/rays/RayTracePropagator.m`
- `src/propagation/rays/OpticalRay.m`
- `src/propagation/rays/CylindricalRay.m`
- `src/propagation/rays/RayBundle.m`
- `src/propagation/rays/RayTracer.m`
