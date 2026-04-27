# Proposal: Octave CI Failure Propagation

## Intent

Ensure Octave CI fails when the portable test runner reports failing tests. Today the workflow pipes `portable_runner()` output to `tee`, but does not assert the returned failure count; this can hide red tests.

## Scope

### In Scope
- Update Octave portable test workflow command to check `portable_runner()` return value.
- Preserve log capture through `tee` and artifact upload.
- Keep MATLAB workflow behavior unchanged because it already checks status.

### Out of Scope
- Refactoring the full MATLAB/Octave test runner.
- Changing test coverage or test list.
- Solving local runtime/license availability.

## Capabilities

### New Capabilities
- `ci-test-status`: CI jobs MUST propagate portable runner failures as non-zero job failures.

### Modified Capabilities
- None.

## Approach

Use the existing MATLAB CI pattern in the Octave workflow: assign `status = portable_runner();` and call `error(...)` when status is non-zero. Keep `set -o pipefail` so the pipeline still fails even with `tee`.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `.github/workflows/octave.yml` | Modified | Portable suite command asserts returned failure count. |
| `tests/portable_runner.m` | Unchanged | Continue returning `totalFailed`; no CLI exit side effect added. |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Quoting mistake in Octave eval command | Medium | Keep command short and mirror MATLAB logic. |
| Existing hidden failures surface in CI | Medium | Desired behavior; fix tests separately if exposed. |

## Rollback Plan

Revert the `.github/workflows/octave.yml` command to the previous direct `portable_runner();` invocation.

## Dependencies

- Working Octave runtime in GitHub Actions or local user environment.

## Success Criteria

- [ ] Octave CI command errors when `portable_runner()` returns nonzero.
- [ ] Logs are still captured in `portable-tests.log`.
- [ ] MATLAB CI remains unchanged.
