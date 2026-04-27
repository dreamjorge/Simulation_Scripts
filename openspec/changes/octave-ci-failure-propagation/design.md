# Design: Octave CI Failure Propagation

## Technical Approach

Keep the runner as a pure function returning `totalFailed`. Fix only the Octave GitHub Actions command so CI interprets that return value as pass/fail, matching the MATLAB workflow.

## Architecture Decisions

### Decision: Assert status in workflow, not runner

**Choice**: Update `.github/workflows/octave.yml` to call `status = portable_runner(); if status ~= 0, error(...); end`.
**Alternatives considered**: Uncomment `exit(1)` in `tests/portable_runner.m`.
**Rationale**: `portable_runner()` is used from MATLAB, Octave, scripts, and interactive sessions. Keeping exit behavior at CI boundary avoids surprising local callers.

### Decision: Preserve log pipeline

**Choice**: Keep `2>&1 | tee portable-tests.log` with `set -o pipefail`.
**Alternatives considered**: Remove `tee` and rely on console logs only.
**Rationale**: The project already uploads logs as artifacts; preserving this keeps diagnostics intact.

## Data Flow

```text
portable_runner() ──returns totalFailed──> Octave eval command
        │                                      │
        └──── stdout/stderr ──────────────────┴──> tee portable-tests.log
                                               │
                                    status ~= 0 -> error -> CI fails
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `.github/workflows/octave.yml` | Modify | Check portable runner status in Octave portable-tests job. |
| `.github/workflows/matlab.yml` | None | Already checks status. |
| `tests/portable_runner.m` | None | Preserve function-style return behavior. |

## Interfaces / Contracts

`portable_runner()` continues to return numeric `totalFailed`; callers decide whether to raise/exit.

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Static | Workflow command checks status | Inspect `.github/workflows/octave.yml`. |
| Integration | Octave portable suite | Run command with working Octave runtime. |
| Regression | MATLAB CI behavior | Confirm workflow unchanged. |

## Migration / Rollout

No migration required.

## Open Questions

- [ ] Can local verification be run under Windows user `uib95096` non-interactively from this session?
