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
| `GaussianBeam(0, params)` | `ExampleRayTracing.m`, `tests/test_RayTracing.m` | `GaussianBeam(w0, lambda)` | Legacy adapter added to GaussianBeam constructor | ✅ Done |
| `GaussianBeam(R, GaussianParameters)` | `examples/MainGauss_refactored.m`, `examples/MainGauss.m`, `examples/MainAnalyticPropagationGauss.m` | `GaussianBeam(w0, lambda)` + `opticalField(X,Y,z)` | Legacy adapter added to GaussianBeam constructor | ✅ Done |
| Legacy aliases `HankeleHermite`, `HankeleLaguerre` | `ParaxialBeams/Hankele*.m` | `HankelHermite`, `HankelLaguerre` | Emit migration warning on use | ✅ Done |
| Dual constructor in `HankelLaguerre` | `ParaxialBeams/HankelLaguerre.m` | Single modern constructor + separate legacy adapter | Legacy paths retained for backward compat | ✅ Done |
| Plot helpers in `Addons/Plots_Functions/*` | Legacy examples | `VisualizationUtils` unified plotting API | Legacy wrappers retained; VisualizationUtils is modern entrypoint | ✅ Defined |
| Split test runners (`test_all` vs `portable_runner`) | `tests/` | One canonical entrypoint | `test_all.m` delegates to `portable_runner()` | ✅ Done |

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

- [x] Compatibility adapter layer created under `legacy/compat/`. (Adapters embedded in beam classes)
- [x] Canonical examples use modern constructors only.
- [x] `tests/test_all.m` and `tests/portable_runner.m` execute the same suite surface.
- [x] Legacy aliases emit migration warning with replacement path.
- [x] Plot API ownership defined (single modern entrypoint + legacy wrappers).
- [ ] Migration notes added for users running historical scripts. (See below)

## Risks and Mitigations

- Risk: Breaking old thesis/research scripts.
  - Mitigation: keep adapters and run legacy compatibility tests separately.
- Risk: API confusion during transition.
  - Mitigation: one canonical contract in README and architecture docs.
- Risk: Duplicate plotting APIs persist.
  - Mitigation: deprecate direct addon plotting in stages, with wrappers.

## Migration Notes for Legacy Script Users

If you are running historical scripts and encounter warnings or deprecation notices:

### HankeleHermite / HankeleLaguerre Warnings
These aliases are deprecated. Replace:
```matlab
% Old (deprecated)
h = HankeleHermite(x, y, params, type);

% New (recommended)
h = HankelHermite(x, y, params, type);
```

### GaussianBeam Legacy Constructors
The `GaussianBeam(R, params)` and `GaussianBeam(X, Y, params)` forms are deprecated:
```matlab
% Old (deprecated)
gb = GaussianBeam(0, gaussianParams);

% New (recommended)
gb = GaussianBeam(w0, lambda);
field = gb.opticalField(X, Y, z);
```

### Plotting
- **Modern entrypoint**: `VisualizationUtils` class (static methods for ray visualization)
- **Legacy wrappers**: `Addons/Plots_Functions/*.m` (retained for backward compatibility)

### Testing Legacy Scripts
Run `tests/run_migration_smoke.m` to verify your legacy scripts still work with the adapter layer.
