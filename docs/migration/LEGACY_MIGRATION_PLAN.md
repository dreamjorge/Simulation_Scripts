# Legacy Migration Plan (Strangler Strategy)

## Goal

Migrate legacy scripts and mixed APIs to a clear, maintainable structure without breaking current research workflows.

This plan is the Week 1 baseline for an incremental migration. It does not redesign physics models; it separates responsibilities and introduces compatibility adapters.

## Scope (Week 1)

- Freeze the modern public API contract.
- Document legacy-to-modern mappings.
- Classify examples into canonical vs legacy.
- Define acceptance gates to start Week 2 compatibility adapters.

## Out of Scope (Week 1)

- Full package migration to `+paraxial/...`.
- Removing legacy files.
- Rewriting all historical examples.
- Deep OO redesign of beam internals.

## Public API Contract (Frozen)

- Beam constructor contract (modern):
  - `GaussianBeam(w0, lambda)`
  - `HermiteBeam(w0, lambda, n, m)`
  - `LaguerreBeam(w0, lambda, l, p)`
- Beam field entrypoint: `opticalField(X, Y, z)`
- Beam metadata:
  - `getParameters(z)`
  - `beamName()`
- Propagation strategies:
  - `FFTPropagator(grid, lambda[, z0])`
  - `AnalyticPropagator(grid[, z0])`
  - `RayTracePropagator(grid[, method, dz, z0])`

## Legacy -> Modern Compatibility Matrix

| Legacy Surface | Current Usage | Modern Replacement | Adapter Action (Week 2) | Status |
|---|---|---|---|---|
| `GaussianBeam(0, params)` | `ExampleRayTracing.m`, `tests/test_RayTracing.m` | `GaussianBeam(w0, lambda)` | Add adapter helper or explicit migration in examples/tests | needs-fix |
| `GaussianBeam(R, GaussianParameters)` | `examples/MainGauss_refactored.m`, `examples/MainGauss.m`, `examples/MainAnalyticPropagationGauss.m` | `GaussianBeam(w0, lambda)` + `opticalField(X,Y,z)` | Keep legacy scripts, migrate canonical scripts first | mixed |
| Legacy aliases `HankeleHermite`, `HankeleLaguerre` | `ParaxialBeams/Hankele*.m` | `HankelHermite`, `HankelLaguerre` | Move aliases under `legacy/compat` with warning | needs-fix |
| Dual constructor in `HankelLaguerre` | `ParaxialBeams/HankelLaguerre.m` | Single modern constructor + separate legacy adapter | Extract legacy constructor path into adapter | needs-fix |
| Plot helpers in `Addons/Plots_Functions/*` | Legacy examples | `VisualizationUtils` unified plotting API | Wrap old plot functions and deprecate direct use | mixed |
| Split test runners (`test_all` vs `portable_runner`) | `tests/` | One canonical entrypoint | Make `test_all.m` delegate to `portable_runner()` | needs-fix |

## Example Classification

### Canonical (entrypoint for new users)

- `examples/MainGauss_refactored.m` (requires API cleanup in Week 2)
- `examples/MainMultiMode.m`
- `ExampleRayTracing.m` (requires API cleanup in Week 2)

### Legacy (kept for reproducibility)

All other scripts in `examples/` are legacy research or historical variants and remain available, but should not be used as onboarding entrypoints.

Known legacy-marked scripts: 32

## Proposed Target Structure (post-migration)

```text
Simulation_Scripts/
  src/
    beams/
    parameters/
    propagation/
      field/
      rays/
    visualization/
  legacy/
    compat/
    examples/
  examples/
    canonical/
  tests/
    modern/
    legacy_compat/
```

## Week 2 Acceptance Checklist

- [ ] Compatibility adapter layer created under `legacy/compat/`.
- [ ] Canonical examples use modern constructors only.
- [ ] `tests/test_all.m` and `tests/portable_runner.m` execute the same suite surface.
- [ ] Legacy aliases emit migration warning with replacement path.
- [ ] Plot API ownership defined (single modern entrypoint + legacy wrappers).
- [ ] Migration notes added for users running historical scripts.

## Risks and Mitigations

- Risk: Breaking old thesis/research scripts.
  - Mitigation: keep adapters and run legacy compatibility tests separately.
- Risk: API confusion during transition.
  - Mitigation: one canonical contract in README and architecture docs.
- Risk: Duplicate plotting APIs persist.
  - Mitigation: deprecate direct addon plotting in stages, with wrappers.
