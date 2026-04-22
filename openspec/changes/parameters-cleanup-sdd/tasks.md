# Tasks: parameters-cleanup-sdd

## Phase 1: Baseline and test scaffolding

- [x] 1.1 Create `tests/modern/test_BeamComputation.m` with failing checks for core formulas
- [x] 1.2 Run targeted test to confirm failure due to missing class
- [x] 1.3 Capture baseline behavior for `tests/modern/test_GaussianParameters.m`

## Phase 2: Computation layer introduction

- [x] 2.1 Create `src/computation/BeamComputation.m` with core static methods
- [x] 2.2 Run `test_BeamComputation.m` and make it pass
- [x] 2.3 Add/adjust tolerances for Octave/MATLAB numeric consistency

## Phase 3: Parameter delegation refactor

- [x] 3.1 Refactor `GaussianParameters` to delegate formulas to `BeamComputation`
- [x] 3.2 Refactor `HermiteParameters` formula methods to delegation
- [x] 3.3 Refactor `LaguerreParameters` formula methods to delegation
- [x] 3.4 Refactor `ElegantHermiteParameters` to delegation
- [x] 3.5 Refactor `ElegantLaguerreParameters` to delegation
- [x] 3.6 Update/add tests in `tests/modern/test_*Parameters.m` for equivalence

## Phase 4: Legacy helper extraction

- [x] 4.1 Create `src/computation/HermiteComputation.m` with migrated implementation
- [x] 4.2 Keep `HermiteParameters.getHermiteSolutions(...)` as compatibility shim
- [x] 4.3 Add compatibility assertion in `tests/legacy_compat/`

## Phase 5: Path and runner integration

- [x] 5.1 Add `src/computation` to `setpaths.m`
- [x] 5.2 Add `src/computation` to `tests/portable_runner.m`
- [x] 5.3 Verify targeted tests can resolve new classes without manual addpath

## Phase 6: Full verification

- [x] 6.1 Run `tests/modern/test_BeamComputation.m`
- [x] 6.2 Run `tests/test_all.m`
- [x] 6.3 Run `tests/legacy_compat/run_legacy_compat.m`
- [x] 6.4 Smoke run canonical examples for API compatibility

## Phase 7: Documentation and closeout

- [x] 7.1 Update `docs/ARCHITECTURE.md` dependency direction
- [x] 7.2 Document compatibility decisions and deferred items
- [x] 7.3 Final review of changed files and unresolved risks
