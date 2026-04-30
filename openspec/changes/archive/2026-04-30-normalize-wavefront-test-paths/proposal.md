# Proposal: Normalize Wavefront Test Paths

## Intent

Remove misleading path warnings from `tests/modern/test_Wavefront.m`. The full suite passes, but this test computes `repoRoot` as `tests/modern/..`, causing failed `addpath()` attempts for non-existent `tests/src/*` and `tests/ParaxialBeams/*` paths. Passing tests with false warnings are noise; noise hides real failures.

## Scope

### In Scope
- Fix `test_Wavefront.m` repo-root/path initialization.
- Preserve portable execution through `tests/portable_runner.m` and direct script execution.
- Add/extend guardrails so modern tests do not introduce stale `tests/modern/..` root assumptions.

### Out of Scope
- Wavefront numerical behavior changes.
- Beam physics, propagation, or Zernike formula changes.
- Package artifact builds.
- Broad `src/` compatibility reduction.

## Capabilities

### New Capabilities
- `wavefront-test-path-policy`: Wavefront tests use repo-root-aware setup without stale path warnings.

### Modified Capabilities
- `runner-path-policy`: Direct modern test invocation SHOULD preserve package-parent semantics consistently with the portable runner.

## Approach

Use the same repo-root semantics already documented for the portable runner: derive the actual repo root from `tests/modern/../..`, add repo root for `+paraxial/`, and keep any compatibility paths explicit. Add a stable guardrail that prevents modern tests from treating `tests/` as the repo root.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `tests/modern/test_Wavefront.m` | Modified | Correct repo root derivation and path setup comments. |
| `tests/modern/test_RepositoryGuardrails.m` | Modified | Add structural check for stale modern-test path roots. |
| `openspec/specs/runner-path-policy/spec.md` | Modified | Extend direct test invocation expectations. |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Direct Wavefront test loses needed paths | Low | Run direct Wavefront test and full portable suite. |
| Guardrail becomes brittle | Low | Check stable path anti-patterns, not prose formatting. |

## Rollback Plan

Revert the test path change and guardrail addition. No runtime code or physics behavior changes are involved.

## Dependencies

- Existing `tests/portable_runner.m` path policy.
- Existing `tests/modern/test_Wavefront.m` script.

## Success Criteria

- [ ] Direct Wavefront test runs without missing-path `addpath` warnings.
- [ ] Full portable Octave suite passes.
- [ ] Guardrails pass and detect stale modern-test repo-root assumptions.
