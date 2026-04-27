# Design: legacy-alias-removal-sdd

## Context

Legacy aliases currently provide constructor/static delegation to modern
`HankelHermite` and `HankelLaguerre`. They are deprecated and documented for
future removal. We need a safe transition that avoids accidental breakage in
research scripts and CI.

## Architecture Decision

Use a **gate-driven deprecation removal**:

1. Verify readiness gates (Usage/Test/Docs/Release) with explicit checks.
2. Apply removal in a single focused commit.
3. Convert alias compatibility tests into migration assertions (or archive
   alias-only suites) in the same change.

This keeps behavior deterministic and rollback simple.

## Implementation Plan

### Phase A — Gate Verification

- Confirm no `Hankele*` usage in `examples/canonical` and `tests/modern`
  (already enforced by guardrail).
- Scan repository for remaining internal references and classify as expected
  (legacy docs/tests) vs unexpected.
- Capture evidence in `verify-report.md`.

### Phase B — Alias Removal Patch

- Delete:
  - `legacy/compat/HankeleHermite.m`
  - `legacy/compat/HankeleLaguerre.m`
- Update path/bootstrap docs if they mention loading those aliases.

### Phase C — Test Strategy Update

Two valid approaches (choose one during implementation):

- **Option 1 (preferred):** Replace alias constructor tests with migration
  behavior tests that assert clear errors/messages and replacement path.
- **Option 2:** Archive alias-only tests and keep modernization guardrails +
  full portable suite as acceptance gates.

### Phase D — Documentation + Release Notes

- Update `LEGACY_MIGRATION_PLAN.md` to reflect completed removal.
- Add release checkpoint with commit hash + rollback command.

## Risks and Mitigations

- **Risk:** hidden external scripts still call `Hankele*`.
  - **Mitigation:** require one release-cycle warning period evidence before
    removal; publish replacement snippet prominently.

- **Risk:** CI false confidence if legacy tests are removed incorrectly.
  - **Mitigation:** keep `portable_runner` + guardrail + targeted migration
    assertion tests as mandatory checks.

## Rollback

Rollback is a single revert of the alias-removal commit, restoring both files
and previous legacy compatibility tests.
