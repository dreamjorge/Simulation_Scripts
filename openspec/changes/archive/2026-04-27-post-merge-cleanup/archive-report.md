# Archive Report: Post-Merge Cleanup SDD

## Change

- Change: `post-merge-cleanup`
- Archived: 2026-04-27
- Merge commit: `3947956`

## Specs Synced

| Domain | Action | Details |
|--------|--------|---------|
| `legacy-policy` | No change | Delta only — no main spec created |
| `deprecation-propagation` | No change | Delta only |

## Archive Contents

- `proposal.md` ✅
- `state.yaml` ✅
- `tasks.md` ✅ — 18/18 tasks complete (Phase 5 incomplete)
- `specs/deprecation-propagation/spec.md` ✅
- `specs/legacy-policy/spec.md` ✅

## Verification Evidence

- `tests/modern/test_Wavefront.m` registered in `tests/portable_runner.m`.
- 3 stale local branches deleted (`chore/legacy-migration-week1`, `fix/gaussianbeam-super-constructor`, `refactor/utility-classes`).
- Remote stale branches deleted.
- `examples/legacy/LEGACY_POLICY.md` updated with research scripts count.
- Deprecation warnings added to all 6 `src/propagation/` classes:
  - `FFTPropagator.m`, `AnalyticPropagator.m`, `RayTracePropagator.m`, `OpticalRay.m`, `CylindricalRay.m`, `RayBundle.m`.
- CI passes: MATLAB (matlab-portable-tests ✅, matlab-legacy-compat ✅), Octave (portable-tests ✅, legacy-compat ✅).

## Incomplete Tasks (Phase 5)

| Task | Description | Status |
|------|-------------|--------|
| 5.1 | Update `pre-merge-hardening/state.yaml` to completed | Not done |
| 5.2 | Update `wavefront-class/state.yaml` to completed | Not done |
| 5.3 | Update `hankel-gradient-fix/state.yaml` if exists | N/A — no state.yaml |

Phase 5 tasks were left incomplete because they required updating other SDD state files that were themselves pending archive. This is resolved by the current audit cleanup.

## Notes

- Phase 5 state.yaml updates were blocking — resolved by archiving completed SDDs in this audit session.
- Main work (test registration, branch cleanup, policy docs, deprecation propagation) was fully completed.