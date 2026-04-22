# Alias Removal Release Plan

## Objective

Define a named, non-ad-hoc release milestone for removing deprecated aliases:

- `legacy/compat/HankeleHermite.m`
- `legacy/compat/HankeleLaguerre.m`

## Milestone

- **Milestone ID:** `legacy-alias-removal-r1`
- **Target Tag:** `v2026.05-legacy-alias-removal`
- **Scope:** remove legacy aliases + transition alias-only legacy tests to migration assertions

## Warning-Cycle Evidence Policy

A "stable release cycle" is satisfied when deprecation warnings are observed and
documented across at least two consecutive release checkpoints while modern and
canonical surfaces remain alias-free.

### Evidence snapshots

1. `docs/migration/RELEASE_CHECKPOINT_2026-04-15.md`
   - deprecation warnings documented
2. `docs/migration/RELEASE_CHECKPOINT_2026-04-22.md`
   - post-merge hardening complete, guardrail enforced
3. `openspec/changes/legacy-alias-removal-sdd/verify-report.md`
   - guardrail PASS, portable PASS, no alias usage in protected modern surfaces

## Breaking-Change Note Template (for release notes)

```markdown
### Breaking Change: Legacy Hankele* aliases removed

The deprecated aliases `HankeleHermite` and `HankeleLaguerre` were removed.

Use:
- `HankelHermite` instead of `HankeleHermite`
- `HankelLaguerre` instead of `HankeleLaguerre`

If you maintain historical scripts, update imports/usages before upgrading to
`v2026.05-legacy-alias-removal`.
```

## Rollback

If emergency rollback is required after release:

```bash
git revert <alias-removal-commit>
```

And restore release notes to indicate rollback status.
