# Parameters Cleanup Architecture Design (Approach A)

Date: 2026-04-21  
Status: Proposed (approved for planning)

## 1) Goal

Reduce fragility and improve extensibility by separating **parameter data modeling** from **optical computation logic**, while preserving:

- Octave/MATLAB compatibility (Octave 11.1.0+, MATLAB R2020b+)
- no external dependencies
- current public API behavior and signatures

## 2) Context and Problem

Current `src/parameters/*` classes mix responsibilities:

- data/snapshot state (`zCoordinate`, `InitialWaist`, `Lambda`, etc.)
- dynamic numerical computations (`waist(z)`, `gouyPhase(z)`, `radius(z)`, `computeAlphaAtZ(z)`)
- legacy/static helpers in parameter classes (e.g. `HermiteParameters.getHermiteSolutions`)

This coupling increases fragility and makes extensions inconsistent across beam families.

## 3) Design Principles

1. **Single Responsibility**: parameter objects represent model/state; computation utilities perform formulas.
2. **Backwards Compatibility First**: keep public API stable for beams, propagators, tests, and legacy scripts.
3. **Incremental Migration**: delegation-first refactor; no big-bang rewrite.
4. **Physics Consistency**: preserve existing formulas and phase conventions.

## 4) Proposed Architecture

### 4.1 New computation layer

Add a new package folder:

- `src/computation/BeamComputation.m`
- `src/computation/HermiteComputation.m` (initially for legacy Hermite helper migration)

`BeamComputation` will provide stateless static functions:

- `rayleighDistance(w0, lambda)`
- `waveNumber(lambda)`
- `waist(w0, z, lambda, zr)`
- `gouyPhase(z, zr)`
- `radiusOfCurvature(z, zr)`
- `complexBeamParameter(z, zr, k)`
- optional higher-order helpers (Hermite/Laguerre waist and modal Gouy)

### 4.2 Parameter classes remain public façade

Keep existing parameter classes and signatures:

- `GaussianParameters`
- `HermiteParameters`
- `LaguerreParameters`
- `ElegantHermiteParameters`
- `ElegantLaguerreParameters`

They remain API façade/backward-compat wrappers and **delegate computations** to `BeamComputation`.

### 4.3 Legacy helper relocation

Move `HermiteParameters.getHermiteSolutions(...)` internals to:

- `HermiteComputation.hermiteSolutions(...)`

Keep compatibility shim in `HermiteParameters`:

- `getHermiteSolutions(...)` delegates and is documented as legacy/deprecated.

## 5) Dependency Direction (Target)

Target flow:

`Beams/Propagators -> Parameters (API façade) -> BeamComputation (pure formulas)`

Compatibility utilities remain available:

`PhysicalConstants` can continue existing; if needed, it may delegate to `BeamComputation` internally without API changes.

## 6) Non-Goals

- No package migration to `+paraxial/...` in this change
- No change to `ParaxialBeam` public contract
- No change to `BeamFactory` signatures
- No rewrite of beam optical field implementations beyond necessary delegation updates

## 7) Migration Strategy

### Phase 1 — Introduce computation utilities

- Add `src/computation/BeamComputation.m`
- Add tests for formulas in isolation

### Phase 2 — Delegation in parameter classes

- Update dynamic/snapshot methods in all parameter classes to call `BeamComputation`
- Keep existing method/property names and behavior

### Phase 3 — Move legacy static helper

- Add `HermiteComputation.m`
- Keep compatibility shim in `HermiteParameters`

### Phase 4 — Verify integration

- Run portable suite and legacy compatibility suite
- Verify no API regressions in canonical examples

## 8) Testing Strategy

### Unit tests (new)

- `BeamComputation` formulas:
  - `w(0) = w0`
  - `w(zR) = w0*sqrt(2)`
  - `R(0) = Inf`
  - `k = 2*pi/lambda`
  - `alpha(z) = i*k/(2*(z+i*zR))`

### Regression tests (existing)

- `tests/test_all.m` (portable wrapper)
- `tests/legacy_compat/run_legacy_compat.m`
- edge-case suites already present under `tests/edge_cases/`

## 9) Risks and Mitigations

### Risk A: Silent behavior drift from formula refactor

Mitigation:

- keep formulas mathematically identical
- add numerical equivalence tests with tolerances

### Risk B: Legacy script breakage for static helpers

Mitigation:

- preserve shim methods in parameter classes
- document deprecation path, do not remove in this phase

### Risk C: Path/addpath issues for new folder

Mitigation:

- update path setup (`setpaths.m` and test runner paths) to include `src/computation`
- verify both Octave and MATLAB

## 10) Acceptance Criteria

1. Public APIs and signatures unchanged for beams/parameters/factory.
2. Parameter classes delegate calculations to computation utilities.
3. `HermiteParameters.getHermiteSolutions` preserved via delegation shim.
4. Portable + legacy suites pass.
5. Documentation updated with the new dependency direction.

## 11) Follow-up Opportunities

After this change stabilizes:

- define stricter input validation boundaries (constructor vs computation layer)
- evaluate gradual namespace migration (`+paraxial/...`)
- unify duplicated beam formula fragments using shared computation helpers
