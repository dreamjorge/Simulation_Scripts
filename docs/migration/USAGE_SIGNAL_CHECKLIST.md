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

- [ ] No unresolved user reports requesting continued support for `Hankele*` aliases
- [ ] Any reported usage has migration path acknowledged with timeline
- [ ] Breaking-change announcement published before removal merge
- [ ] Observation window completed

## Decision Log

| Date | Signal Source | Finding | Action |
|------|----------------|---------|--------|
| 2026-04-22 | Baseline setup | Checklist initialized | Start observation window |

## Final Sign-off

- [ ] Usage gate approved by maintainer
- Sign-off date: ____
- Signed by: ____
