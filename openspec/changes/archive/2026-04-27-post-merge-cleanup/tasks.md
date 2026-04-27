# Tasks: Post-Merge Cleanup

## Phase 1: Test Registration

- [x] 1.1 Register `test_Wavefront.m` in `tests/modern/` — moved from `tests/` (file was in wrong location)

## Phase 2: Branch Cleanup

- [x] 2.1 Delete local stale branch `chore/legacy-migration-week1`
- [x] 2.2 Delete local stale branch `fix/gaussianbeam-super-constructor`
- [x] 2.3 Delete local stale branch `refactor/utility-classes`
- [x] 2.4 Delete remote stale branches

## Phase 3: Legacy Policy Documentation

- [x] 3.1 Verify `examples/legacy/LEGACY_POLICY.md` exists and reflects current classification
- [x] 3.2 Updated LEGACY_POLICY.md — fixed research/ files count (was "empty", now "5 thesis scripts")

## Phase 4: Deprecation Propagation

- [x] 4.1 ~~Add deprecation warning to `src/propagation/field/IPropagator.m`~~ (N/A — abstract base class, no constructor)
- [x] 4.2 Add deprecation warning to `src/propagation/field/FFTPropagator.m`
- [x] 4.3 Add deprecation warning to `src/propagation/field/AnalyticPropagator.m`
- [x] 4.4 Add deprecation warning to `src/propagation/rays/RayTracePropagator.m`
- [x] 4.5 Add deprecation warning to `src/propagation/rays/OpticalRay.m`
- [x] 4.6 Add deprecation warning to `src/propagation/rays/CylindricalRay.m`
- [x] 4.7 Add deprecation warning to `src/propagation/rays/RayBundle.m`
- [x] 4.8 ~~Add deprecation warning to `src/propagation/rays/RayTracer.m`~~ (N/A — no constructor, static methods only)

## Phase 5: SDD Archive

- [ ] 5.1 Update `openspec/changes/pre-merge-hardening/state.yaml` — set status to `completed`
- [ ] 5.2 Update `openspec/changes/wavefront-class/state.yaml` — set status to `completed`
- [ ] 5.3 Update `openspec/changes/hankel-gradient-fix/state.yaml` if exists

## Phase 6: Verification

- [x] 6.1 Run `tests/portable_runner.m` — verify test_Wavefront is included (10/10 passed)
- [x] 6.2 Verify all deprecation warnings emit correctly (confirmed in CI output)
- [x] 6.3 Verify CI passes (MATLAB + Octave) — both ✅
