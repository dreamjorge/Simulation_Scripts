# Proposal: legacy-alias-removal-sdd

## Why

The repository still ships `legacy/compat/HankeleHermite.m` and
`legacy/compat/HankeleLaguerre.m` as compatibility aliases. This is useful for
historical scripts, but it keeps technical debt alive and blurs the modern API
surface.

After PR #26, we already have:

- readiness gates documented in `docs/migration/LEGACY_MIGRATION_PLAN.md`
- a guardrail test (`test_LegacyAliasGuardrail.m`) that prevents deprecated
  aliases from re-entering canonical/modern surfaces

This change defines a controlled, evidence-based path to remove legacy aliases
without breaking active workflows.

## What Changes

1. Add objective readiness verification tasks for Usage/Test/Docs/Release gates.
2. Prepare alias-removal patch (delete `legacy/compat/Hankele*.m`) behind a
   gate pass decision.
3. Update legacy compatibility tests so they assert migration behavior after
   removal (or archive alias-only tests, depending on chosen strategy).
4. Update migration docs and release checkpoint with exact removal commit and
   rollback instructions.

## Non-Goals

- No migration to `+paraxial` packaging in this track.
- No redesign of beam physics internals.
- No API changes to modern constructors/`opticalField` signatures.

## Success Criteria

- All readiness gates pass with evidence.
- `portable_runner` remains green after alias removal.
- Canonical examples and modern tests remain alias-free.
- Migration docs include explicit "removed in release" notice and replacement
  snippets.
