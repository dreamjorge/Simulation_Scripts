# Release Checkpoint - 2026-04-22

## Scope

This checkpoint closes the post-merge architecture hardening after:

- PR #24 (`feat/parameters-cleanup-sdd`)
- PR #25 (`chore/parameters-next-steps`)

Focus is architecture consistency, migration readiness, and guardrails against legacy reintroduction.

## Included

- `BeamComputation` + `HermiteComputation` extraction is merged in `master`.
- `GaussianParameters` now delegates static helpers through `BeamComputation` (no residual formula split).
- Canonical smoke compatibility validated and plotting issues fixed in `MainGauss_refactored`.
- Legacy removal policy documented with explicit readiness gates (Usage/Test/Docs/Release).
- New guardrail test added:
  - `tests/modern/test_LegacyAliasGuardrail.m`
  - Enforced from `portable_runner.m`
  - Fails if deprecated `Hankele*` aliases appear in:
    - `examples/canonical/`
    - `tests/modern/`

## Validation Snapshot

- `master` includes merge commit `e44d32e` (PR #25).
- `portable_runner` now includes the legacy alias guardrail test.

## Not Included

- Retirement of `tests/legacy_compat/`.
- Package migration from `src/` to `+paraxial/`.

## Next Phase Recommendation

1. Keep guardrail green for at least one release cycle after alias removal.
2. Monitor post-removal support signals via `USAGE_SIGNAL_CHECKLIST.md`.
3. Maintain rollback readiness for removal commit.

## Alias Removal Execution (Completed)

- Removed files:
  - `legacy/compat/HankeleHermite.m`
  - `legacy/compat/HankeleLaguerre.m`
- Removal commit hash: `58b02b0`

Rollback:

```bash
git revert 58b02b0
```
