# Archive Report: Legacy Alias Removal SDD

## Change

- Change: `legacy-alias-removal-sdd`
- Archived: 2026-04-27
- Merge commit: `58b02b0` (direct to master)

## Specs Synced

| Domain | Action | Details |
|--------|--------|---------|
| `legacy-migration-gates` | No change | Delta only — no main spec created |
| `legacy-alias-removal` | No change | Delta only |

## Archive Contents

- `proposal.md` ✅
- `design.md` ✅
- `state.yaml` ✅
- `tasks.md` ✅ — all phases complete
- `verify-report.md` ✅ — verdict PASS
- `specs/legacy-alias-removal/spec.md` ✅
- `specs/legacy-migration-gates/spec.md` ✅

## Verification Evidence

- PR #31/34 combined: `HankeleHermite.m` and `HankeleLaguerre.m` removed from `legacy/compat/`.
- `tests/modern/test_LegacyAliasGuardrail.m` passes (no deprecated aliases in modern code).
- `tests/portable_runner.m` green (31/0) after removal.
- Temporary alias-removal test run confirmed full suite passes without aliases.
- `LEGACY_ALIAS_REMOVAL_MODE=1` flag enables legacy tests to pass post-removal.

## Source of Truth Updated

- `legacy/compat/HankeleHermite.m` — removed
- `legacy/compat/HankeleLaguerre.m` — removed
- `docs/migration/LEGACY_MIGRATION_PLAN.md` — updated with removal status

## Notes

- Alias removal was executed before Phase 4 archive tasks were completed; archive created retroactively.
- Migration gate specs preserved for future reference.