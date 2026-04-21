# Proposal: parameters-cleanup-sdd

## Intent

Implement the approved post-merge architecture direction to reduce fragility and improve extensibility by separating parameter data modeling from numerical beam computations, while preserving MATLAB/Octave compatibility and the current public API.

## Scope

### In Scope
- Introduce a stateless computation layer under `src/computation/`
- Delegate formulas currently embedded in `src/parameters/*.m` to computation utilities
- Preserve current constructor/method/property signatures in parameter classes
- Extract legacy Hermite helper logic to dedicated computation utility with compatibility shim
- Update path setup and test runner to include `src/computation`
- Add/adjust tests for formula equivalence and backward compatibility

### Out of Scope
- Namespace/package migration to `+paraxial/...`
- Public API redesign for beams/propagators/factory
- Major refactor of beam field equations beyond delegation needs
- Removal of legacy compatibility entry points

## Capabilities

### New Capabilities
- `parameter-computation-layer`: Centralized stateless formulas for Gaussian/Hermite/Laguerre/elegant parameter calculations

### Modified Capabilities
- `parameter-model-delegation`: Existing parameter classes become API façades that delegate formula evaluation
- `legacy-hermite-helper-compat`: Preserve `HermiteParameters.getHermiteSolutions(...)` via compatibility shim

## Approach

1. Add `BeamComputation` with pure static math functions.
2. Refactor parameter classes to delegate dynamic and snapshot calculations.
3. Extract Hermite legacy helper to `HermiteComputation` and keep backward-compatible delegation shim.
4. Wire `src/computation` into path setup (`setpaths.m`, `tests/portable_runner.m`).
5. Validate with focused formula tests plus full portable and legacy suites.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `src/computation/` | New | Stateless computation utilities |
| `src/parameters/*.m` | Modified | Delegation to computation layer, API preserved |
| `setpaths.m` | Modified | Add computation folder to runtime path |
| `tests/portable_runner.m` | Modified | Add computation folder to test path |
| `tests/modern/*` | Modified/New | Formula equivalence and delegation tests |
| `tests/legacy_compat/*` | Modified | Compatibility validation for Hermite helper |
| `docs/ARCHITECTURE.md` | Modified | Dependency direction update |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Silent numerical drift in formulas | Medium | Golden formula tests with tolerances |
| Legacy script breakage | Medium | Keep shim methods and run legacy suite |
| Path-resolution failures in Octave/MATLAB | Medium | Update both setpaths and portable runner |

## Success Criteria

- [ ] `src/computation/BeamComputation.m` added and covered by tests
- [ ] Parameter classes delegate formulas without API/signature changes
- [ ] Legacy Hermite helper remains functional through shim
- [ ] Portable suite passes (`tests/test_all.m`)
- [ ] Legacy compatibility suite passes (`tests/legacy_compat/run_legacy_compat.m`)
- [ ] Architecture docs updated with new dependency direction
