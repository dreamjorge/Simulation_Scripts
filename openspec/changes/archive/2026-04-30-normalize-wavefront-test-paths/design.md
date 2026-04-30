# Design: Normalize Wavefront Test Paths

## Technical Approach

Make `test_Wavefront.m` follow the existing package-parent path policy. The script currently assumes `fullfile(testDir, '..')` is the repo root, but from `tests/modern/` that resolves to `tests/`. The fix is to derive `repoRoot = fullfile(testDir, '..', '..')`, add repo root for `+paraxial/`, then add compatibility/utility paths from that root.

## Architecture Decisions

| Decision | Choice | Alternatives considered | Rationale |
|----------|--------|-------------------------|-----------|
| Repo root derivation | Use `tests/modern/../..` | Use current `..` | Current code points at `tests/`, producing false warnings. |
| Package semantics | Add repo root | Add internal `+paraxial` folders | Matches MATLAB/Octave package semantics and runner policy. |
| Guardrail | Scan for stale modern-test path setup | Rely on warnings/manual review | Warnings are easy to ignore; guardrails keep the invariant executable. |

## Data Flow

```text
test_Wavefront.m
  └─ mfilename('fullpath')
      └─ tests/modern
          └─ ../.. -> repo root
              ├─ +paraxial package parent
              ├─ src/* deprecated compatibility paths
              └─ ParaxialBeams utilities
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `tests/modern/test_Wavefront.m` | Modify | Correct repo root and path setup. |
| `tests/modern/test_RepositoryGuardrails.m` | Modify | Add anti-pattern check for stale modern-test root assumptions. |
| `openspec/changes/normalize-wavefront-test-paths/specs/*/spec.md` | Create | Capture path behavior requirements. |

## Interfaces / Contracts

No public MATLAB API changes. This affects test setup only.

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Direct script | `test_Wavefront.m` has no missing-path warnings | Run with Octave CLI directly. |
| Guardrail | Stale modern-test root assumptions are detected | Run `test_RepositoryGuardrails.m`. |
| Integration | Full suite remains green | Run `tests/test_all.m`. |

## Migration / Rollout

No migration required.

## Open Questions

None.
