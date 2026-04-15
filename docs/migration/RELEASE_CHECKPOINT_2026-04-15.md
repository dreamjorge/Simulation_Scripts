# Release Checkpoint - 2026-04-15

## Scope

This checkpoint closes the current legacy migration hardening cycle based on the Strangler strategy.
It focuses on backward compatibility and testability, not package migration to `+paraxial`.

## Included

- Legacy alias compatibility retained via `legacy/compat`.
- Static propagation API delegation validated for:
  - `HankeleHermite.getPropagateCartesianRays`
  - `HankeleLaguerre.getPropagateCylindricalRays`
- Edge-case delegation parity tests added (mixed slopes, negative radial scenarios).
- Migration docs and test docs updated to reflect current guarantees.
- Legacy-only compatibility runner added for faster feedback:
  - `run('tests/legacy_compat/run_legacy_compat.m')`

## Validation Snapshot

- `portable_runner()` passes with 0 failures.
- `run_legacy_compat.m` passes all legacy suites.

## Not Included

- Migration of classes from `src/` to `+paraxial/`.
- Removal of legacy aliases.
- Redesign of physics internals.

## Next Phase Recommendation

1. Add CI jobs for:
   - full portable suite (`portable_runner`)
   - legacy-only suite (`run_legacy_compat`)
2. Start a controlled RFC/spec for `+paraxial` packaging with adapter-first rollout.

## Branch Protection Recommendation

For `master`, enable branch protection with:

- Require pull request before merge.
- Require status checks to pass:
  - `portable-tests`
  - `legacy-compat`
- Require branches to be up to date before merging.
- Restrict direct pushes to `master`.
