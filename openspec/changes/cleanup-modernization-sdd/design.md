# Design: Cleanup & Modernization Alignment

## Technical Approach

Perform a bounded cleanup pass: normalize repo/tooling metadata, align docs with current `+paraxial/` and GitHub Actions reality, then add structural guardrail tests. No numerical code or beam formulas change.

## Architecture Decisions

| Decision | Choice | Alternatives considered | Rationale |
|----------|--------|-------------------------|-----------|
| CI ownership | GitHub Actions is canonical; CircleCI removed or marked inactive | Keep both | `.circleci/config.yml` only runs `welcome/run`; dual CI narratives create noise. |
| Tooling metadata | Ignore `.opencode/` unless project intentionally depends on it | Commit package metadata | Current `.opencode/` is untracked local tooling; default should avoid accidental repo noise. |
| Active roadmap | Create/update `docs/ROADMAP.md`; archive root `plan.md` if historical | Keep root plan as active | `plan.md` describes pre-merge hardening already superseded by v2 integration. |
| Guardrails | Add structural tests under `tests/modern/` | Manual review only | Cheap tests prevent docs/examples from drifting back to deprecated surfaces. |

## Data Flow

```text
Developer change
  ├─ updates docs/tooling policy
  ├─ runs portable tests
  └─ guardrail tests inspect docs/examples/factory
        ├─ README / ARCHITECTURE / tests README
        ├─ examples/canonical/*.m
        └─ ParaxialBeams/BeamFactory.m
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `.gitignore` | Modify/Create | Add `.opencode/` if local-tooling policy is ignore. |
| `.circleci/config.yml` | Delete/Archive | Remove stale welcome-only CI config. |
| `README.md` | Modify | Align support matrix, Quick Start, canonical API, release notes. |
| `docs/ARCHITECTURE.md` | Modify | Reflect current `+paraxial/` canonical architecture and deprecated surfaces. |
| `tests/README.md` | Modify | Match actual portable runner commands and CI behavior. |
| `.atl/skill-registry.md` | Modify | Align Octave baseline and test command with README/CI. |
| `docs/ROADMAP.md` | Create | Post-v2 cleanup/modernization roadmap. |
| `plan.md` | Move/Modify | Archive or mark historical. |
| `tests/modern/test_RepositoryGuardrails.m` | Create | Check docs/examples/factory structural invariants. |
| `tests/test_all.m` or `tests/portable_runner.m` | Modify | Register guardrail test. |

## Interfaces / Contracts

No public MATLAB API changes. Guardrail tests are internal validation only.

Stable invariants:
- `+paraxial/` is canonical namespace.
- `BeamFactory.create()` is preferred high-level construction API.
- `src/` is deprecated/transitional.
- GitHub Actions is canonical CI.

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Structural | Canonical examples avoid direct `src/beams` usage | MATLAB/Octave text scan in `test_RepositoryGuardrails.m`. |
| Documentation | README and registry agree on support baseline/API claims | Guardrail assertions on required phrases and forbidden stale claims. |
| Integration | New guardrail test runs in portable suite | Register in existing runner and execute via `tests/test_all.m`. |

## Migration / Rollout

No runtime migration required. Roll out as one cleanup PR. If guardrail tests are too brittle, relax assertions to check only stable invariants.

## Open Questions

- [ ] Should `.opencode/` be ignored or intentionally committed?
- [ ] Is the official Octave baseline 11.1.0+ or lower?
