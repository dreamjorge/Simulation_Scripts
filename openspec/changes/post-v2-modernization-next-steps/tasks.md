# Tasks: Post-v2 Modernization Next Steps

## Phase 1: Runner Path Policy

- [x] 1.1 Update `tests/portable_runner.m` to add repo root under a canonical `+paraxial/` package section.
- [x] 1.2 Relabel `src/*` path additions in `tests/portable_runner.m` as deprecated/transitional compatibility paths.
- [x] 1.3 Confirm `ParaxialBeams/` and `ParaxialBeams/Addons/` remain utility/addon sections, not canonical package paths.

## Phase 2: Guardrails

- [x] 2.1 Extend `tests/modern/test_RepositoryGuardrails.m` to assert `portable_runner.m` does not call `src/*` “modern library paths”.
- [x] 2.2 Add a guardrail that `portable_runner.m` or `setpaths.m` preserves repo-root package parent semantics for `+paraxial/`.
- [x] 2.3 Add roadmap guardrail checks for `docs/ROADMAP.md` and historical root `plan.md` wording.

## Phase 3: Addons Inventory

- [x] 3.1 Create `docs/ADDONS_INVENTORY.md` with classifications for top-level `ParaxialBeams/Addons/*.m` files.
- [x] 3.2 Classify `ParaxialBeams/Addons/Plots_Functions/` as a subdirectory entry with evidence and follow-up action.
- [x] 3.3 Mark uncertain entries as `needs-investigation`; do not delete addon files in this change.

## Phase 4: Compatibility Reduction Planning

- [x] 4.1 Create `docs/COMPATIBILITY_REDUCTION.md` defining gates before reducing `src/` usage or availability.
- [x] 4.2 Update `docs/ROADMAP.md` to reference runner cleanup, addons inventory, and compatibility reduction as active next steps.
- [x] 4.3 Update `tests/README.md` if it conflicts with `setpaths()` or portable runner path policy.

## Phase 5: Verification

- [x] 5.1 Run `octave --no-gui --eval "run('tests/test_all.m')"` and record result. Passed with 33 passed, 0 failed using Octave 11.1.0 absolute path.
- [x] 5.2 Run `octave-cli --no-gui --eval "run('tests/modern/test_RepositoryGuardrails.m')"` and record result. Passed with 12/12 guardrails using Octave 11.1.0 absolute path.
- [x] 5.3 Confirm `git diff` shows no beam physics changes under `+paraxial/+beams/`, `src/beams/`, `src/parameters/`, or `src/computation/`.
