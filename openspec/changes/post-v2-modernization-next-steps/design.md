# Design: Post-v2 Modernization Next Steps

## Technical Approach

Implement a second cleanup wave that strengthens executable policy before deleting anything. The runner keeps compatibility paths but names them correctly; guardrails validate stable invariants; docs capture addons and future `src/` reduction decisions.

## Architecture Decisions

| Decision | Choice | Alternatives considered | Rationale |
|----------|--------|-------------------------|-----------|
| Runner paths | Add repo root for `+paraxial/`; keep `src/*` as deprecated compatibility | Remove `src/*` immediately | Existing legacy compatibility tests still rely on transitional surfaces. |
| Guardrails | Extend `test_RepositoryGuardrails.m` with stable text/path checks | Manual review only | The repo already uses structural guardrails for modernization drift. |
| Addons | Inventory first, no deletion | Delete apparently unused files | Addons may support legacy research/plotting scripts not covered by tests. |
| Compatibility reduction | Create separate plan document | Fold into cleanup PR | Removing `src/` is a public compatibility decision, not hygiene. |

## Data Flow

```text
tests/test_all.m
  └─ portable_runner()
       ├─ repo root              -> +paraxial package resolution
       ├─ ParaxialBeams/Addons   -> utilities + legacy helpers
       ├─ src/*                  -> deprecated compatibility adapters
       └─ tests/*                -> modern, legacy_compat, edge cases

guardrail test
  ├─ reads docs and runner text
  ├─ scans canonical examples
  └─ checks BeamFactory registry
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `tests/portable_runner.m` | Modify | Add canonical package parent path section; relabel `src/*` as deprecated compatibility. |
| `tests/modern/test_RepositoryGuardrails.m` | Modify | Assert runner path policy, roadmap ownership, and future `src/` reduction boundaries. |
| `docs/ROADMAP.md` | Modify | Add active next-wave section and point to this OpenSpec change. |
| `docs/ADDONS_INVENTORY.md` | Create | Classify addon files/directories with evidence and follow-up action. |
| `docs/COMPATIBILITY_REDUCTION.md` | Create | Define migration gates before reducing `src/` dependency. |
| `tests/README.md` | Modify if needed | Prefer `setpaths()` and canonical portable runner commands. |

## Interfaces / Contracts

No public MATLAB API changes. `BeamFactory.supportedTypes()` remains unchanged. `src/` remains available during this change.

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Structural | Runner path labels and package parent path | Text assertions in `test_RepositoryGuardrails.m`. |
| Documentation | Roadmap, addons inventory, compatibility plan exist and use scoped language | Guardrail/manual doc audit. |
| Integration | Cleanup does not break suite | Run `octave --no-gui --eval "run('tests/test_all.m')"`. |

## Migration / Rollout

No runtime migration. This change prepares future migration by documenting gates before `src/` reduction.

## Open Questions

- [ ] Should addon inventory classify nested plotting functions individually or by subdirectory first?
