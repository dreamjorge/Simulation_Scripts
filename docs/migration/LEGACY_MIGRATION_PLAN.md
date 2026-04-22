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
| `GaussianBeam(0, params)` | `examples/canonical/ExampleRayTracing.m`, `tests/modern/test_RayTracing.m` | `GaussianBeam(w0, lambda)` | Legacy adapter added to GaussianBeam constructor | ✅ Done |
| `GaussianBeam(R, GaussianParameters)` | `examples/canonical/MainGauss_refactored.m`, `examples/MainGauss.m`, `examples/MainAnalyticPropagationGauss.m` | `GaussianBeam(w0, lambda)` + `opticalField(X,Y,z)` | Legacy adapter added to GaussianBeam constructor | ✅ Done |
| Legacy aliases `HankeleHermite`, `HankeleLaguerre` | `ParaxialBeams/Hankele*.m` | `HankelHermite`, `HankelLaguerre` | Emit migration warning on use | ✅ Done |
| Legacy alias static propagation API | Historical scripts using `Hankele*.getPropagate*` | `Hankel*.getPropagate*` | Regression tests assert exact delegation parity (normal + edge cases) | ✅ Done |
| Dual constructor in `HankelLaguerre` | `ParaxialBeams/HankelLaguerre.m` | Single modern constructor + separate legacy adapter | Legacy paths retained for backward compat | ✅ Done |
| Plot helpers in `Addons/Plots_Functions/*` | Legacy examples | `VisualizationUtils` unified plotting API | Legacy wrappers retained; VisualizationUtils is modern entrypoint | ✅ Defined |
| Split test runners (`test_all` vs `portable_runner`) | `tests/` | One canonical entrypoint | `test_all.m` delegates to `portable_runner()` | ✅ Done |

## Example Classification

### Canonical (entrypoint for new users)

- `examples/canonical/MainGauss_refactored.m` - uses modern API ✅
- `examples/canonical/MainMultiMode.m` - uses modern API ✅
- `examples/canonical/ExampleRayTracing.m` - migrated to modern API ✅

### Legacy (kept for reproducibility)

All other scripts in `examples/` are legacy research or historical variants and remain available, but should not be used as onboarding entrypoints.

Known legacy-marked scripts: 32

## Current Structure (post-Week 7)

```text
Simulation_Scripts/
  src/
    beams/           % ✅ Beam classes (ParaxialBeam, GaussianBeam, etc.)
    parameters/      % ✅ Parameter classes (GaussianParameters, etc.)
    propagation/
      field/         % ✅ Field propagators (FFT, Analytic)
      rays/          % ✅ Ray propagation (RayTrace, RayBundle, etc.)
    visualization/   % ✅ VisualizationUtils
  ParaxialBeams/    % Utilities (PhysicalConstants, GridUtils, FFTUtils, etc.)
  legacy/
    compat/          % ✅ Legacy aliases (HankeleHermite, HankeleLaguerre)
  examples/
    canonical/       % ✅ Canonical examples for new users
  tests/
    modern/          % ✅ Modern API tests
    legacy_compat/   % ✅ Legacy compatibility tests
  +paraxial/         % 🔜 Future: MATLAB package structure (placeholder)
```

## Week 2 Acceptance Checklist

- [x] Compatibility adapter layer created under `legacy/compat/`. (Adapters embedded in beam classes)
- [x] Canonical examples use modern constructors only.
- [x] `tests/test_all.m` and `tests/portable_runner.m` execute the same suite surface.
- [x] Legacy aliases emit migration warning with replacement path.
- [x] Plot API ownership defined (single modern entrypoint + legacy wrappers).
- [x] Migration notes added for users running historical scripts. (See below)

## Week 3: Implement src/ Structure

- [x] Create `src/` directory structure.
- [x] Move beam classes to `src/beams/`.
- [x] Move parameter classes to `src/parameters/`.
- [x] Move field propagators to `src/propagation/field/`.
- [x] Move ray propagation to `src/propagation/rays/`.
- [x] Move visualization to `src/visualization/`.
- [x] Create `setpaths.m` for easy path configuration.
- [x] Update canonical examples for new structure.
- [x] Update test runners for new structure.

## Week 4: Create legacy/compat Layer

- [x] Create `legacy/compat/` directory structure.
- [x] Move legacy aliases (HankeleHermite, HankeleLaguerre) to `legacy/compat/`.
- [x] Add `legacy/compat/README.md` with usage instructions.
- [x] Update addpath references in legacy examples that use Hankele*.
- [x] Add static alias delegation regression tests (`test_HankelAliasStaticDelegation`, `test_HankelAliasEdgeCases`).

## Week 5: Canonical Examples Folder

- [x] Create `examples/canonical/` directory.
- [x] Move `MainGauss_refactored.m` to `examples/canonical/`.
- [x] Move `MainMultiMode.m` to `examples/canonical/`.
- [x] Move `ExampleRayTracing.m` to `examples/canonical/`.
- [x] Update README with new examples path.

## Week 6: Test Folder Reorganization

- [x] Create `tests/modern/` directory.
- [x] Move modern tests to `tests/modern/`.
- [x] Create `tests/legacy_compat/` directory.
- [x] Move `test_HankelCompatibility.m` to `tests/legacy_compat/`.
- [x] Update `portable_runner.m` for new test structure.

## Week 7: Package Migration (Future Work)

- [x] Create `+paraxial/` package directory structure.
- [x] Add `+paraxial/README.md` with migration notes and planned structure.
- [ ] **Future**: Migrate classes from `src/` to `+paraxial/` package.
- [ ] **Future**: Update imports in examples and tests.
- [ ] **Future**: Deprecate `src/` once `+paraxial/` is stable.

**Note**: Package migration is explicitly **out of scope** for current migration work.
Full package migration will be done in a future phase after `src/` is stable.

## Risks and Mitigations

- Risk: Breaking old thesis/research scripts.
  - Mitigation: keep adapters and run legacy compatibility tests separately.
- Risk: API confusion during transition.
  - Mitigation: one canonical contract in README and architecture docs.
- Risk: Duplicate plotting APIs persist.
  - Mitigation: deprecate direct addon plotting in stages, with wrappers.

## Legacy Removal Readiness Gates (when we can safely remove legacy)

Do **not** remove `legacy/compat/*` or legacy constructor paths until all gates below are green.

### Gate A — Usage Gate

- [x] No internal references to `HankeleHermite` / `HankeleLaguerre` in `src/`, `examples/canonical/`, and `tests/modern/`.
- [ ] At least one migration checkpoint confirms no active external dependency reports from users.

Usage signal execution artifact:

- `docs/migration/USAGE_SIGNAL_CHECKLIST.md`

### Gate B — Test Gate

- [x] `tests/portable_runner.m` passes with legacy aliases still present.
- [x] A dedicated branch run validates portable suite still passes **after** temporarily removing `legacy/compat/Hankele*.m`.
- [x] Legacy-only test suite is either retired or replaced with explicit migration assertions.

### Gate C — Documentation Gate

- [x] README and migration docs no longer recommend legacy aliases anywhere.
- [x] Canonical examples exclusively use modern APIs (`HankelHermite`, `HankelLaguerre`, `opticalField`).
- [x] Release notes include a breaking-change notice and replacement snippets.

Announcement artifact:

- `docs/migration/ALIAS_REMOVAL_ANNOUNCEMENT_TEMPLATE.md`

### Gate D — Release Gate

- [x] Legacy aliases have emitted deprecation warnings for at least one stable release cycle.
- [x] Removal is scheduled in a named release milestone (version/tag), not ad-hoc.

### Removal Checklist (execute only after A+B+C+D)

- [ ] Remove `legacy/compat/HankeleHermite.m` and `legacy/compat/HankeleLaguerre.m`.
- [ ] Remove/update legacy alias tests in `tests/legacy_compat/`.
- [ ] Run full verification (`tests/portable_runner.m` + targeted canonical smokes).
- [ ] Update `docs/migration/RELEASE_CHECKPOINT_*.md` with the exact removal commit hash.

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
