# Tasks: parameters-cleanup-sdd

## Phase 1: Baseline and test scaffolding

- [ ] 1.1 Create `tests/modern/test_BeamComputation.m` with failing checks for core formulas
- [ ] 1.2 Run targeted test to confirm failure due to missing class
- [ ] 1.3 Capture baseline behavior for `tests/modern/test_GaussianParameters.m`

## Phase 2: Computation layer introduction

- [ ] 2.1 Create `src/computation/BeamComputation.m` with core static methods
- [ ] 2.2 Run `test_BeamComputation.m` and make it pass
- [ ] 2.3 Add/adjust tolerances for Octave/MATLAB numeric consistency

## Phase 3: Parameter delegation refactor

- [ ] 3.1 Refactor `GaussianParameters` to delegate formulas to `BeamComputation`
- [ ] 3.2 Refactor `HermiteParameters` formula methods to delegation
- [ ] 3.3 Refactor `LaguerreParameters` formula methods to delegation
- [ ] 3.4 Refactor `ElegantHermiteParameters` to delegation
- [ ] 3.5 Refactor `ElegantLaguerreParameters` to delegation
- [ ] 3.6 Update/add tests in `tests/modern/test_*Parameters.m` for equivalence

## Phase 4: Legacy helper extraction

- [ ] 4.1 Create `src/computation/HermiteComputation.m` with migrated implementation
- [ ] 4.2 Keep `HermiteParameters.getHermiteSolutions(...)` as compatibility shim
- [ ] 4.3 Add compatibility assertion in `tests/legacy_compat/`

## Phase 5: Path and runner integration

- [ ] 5.1 Add `src/computation` to `setpaths.m`
- [ ] 5.2 Add `src/computation` to `tests/portable_runner.m`
- [ ] 5.3 Verify targeted tests can resolve new classes without manual addpath

## Phase 6: Full verification

- [ ] 6.1 Run `tests/modern/test_BeamComputation.m`
- [ ] 6.2 Run `tests/test_all.m`
- [ ] 6.3 Run `tests/legacy_compat/run_legacy_compat.m`
- [ ] 6.4 Smoke run canonical examples for API compatibility

## Phase 7: Documentation and closeout

- [ ] 7.1 Update `docs/ARCHITECTURE.md` dependency direction
- [ ] 7.2 Document compatibility decisions and deferred items
- [ ] 7.3 Final review of changed files and unresolved risks
