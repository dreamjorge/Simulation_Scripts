# Usage Signal Checklist (Legacy Alias Removal)

## Goal

Collect explicit evidence that removing `HankeleHermite` / `HankeleLaguerre`
will not break active users beyond documented migration scope.

## Window

- Start date: 2026-04-22
- Minimum observation window: 14 days (or one stable release cycle)

## Evidence Sources

1. GitHub Issues / Discussions / PR comments
2. Release notes feedback (post-announcement)
3. Internal user reports (if applicable)
4. Direct maintainer contact channel (email/chat) summary

## Commands & Queries

Use these queries to check public signals:

```bash
gh issue list --state all --search "HankeleHermite OR HankeleLaguerre OR Hankel alias"
gh pr list --state all --search "HankeleHermite OR HankeleLaguerre OR Hankel alias"
gh api repos/dreamjorge/Simulation_Scripts/discussions
```

## Acceptance Criteria

- [x] No unresolved user reports requesting continued support for `Hankele*` aliases
- [x] Any reported usage has migration path acknowledged with timeline
- [ ] Breaking-change announcement published before removal merge
- [ ] Observation window completed

## Current Public Signal Snapshot (2026-04-22)

- `gh issue list --state all --search "HankeleHermite OR HankeleLaguerre OR Hankel alias"`
  - Result: no matching issues
- `gh pr list --state all --search "HankeleHermite OR HankeleLaguerre OR Hankel alias"`
  - Result: matched migration PRs only (`#17`, `#18`, `#27`), no unresolved support requests
- `gh api repos/dreamjorge/Simulation_Scripts/discussions`
  - Result: HTTP 410 (Discussions disabled for repo)

## Decision Log

| Date | Signal Source | Finding | Action |
|------|----------------|---------|--------|
| 2026-04-22 | Baseline setup | Checklist initialized | Start observation window |
| 2026-04-22 | GitHub Issues | No alias support requests found | Mark unresolved-request criterion as satisfied |
| 2026-04-22 | GitHub PRs | Only migration-related merged PRs found (#17, #18, #27) | Keep monitoring until window closes |
| 2026-04-22 | GitHub Discussions API | Discussions disabled (HTTP 410) | Treat as N/A source; rely on Issues/PRs |

## Final Sign-off

- [ ] Usage gate approved by maintainer
- Sign-off date: ____
- Signed by: ____
