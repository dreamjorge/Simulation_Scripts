# Tasks: Cleanup & Modernization Alignment

## Phase 1: Repository Hygiene

- [x] 1.1 Decide `.opencode/` policy; update `.gitignore` or commit `.opencode/package.json` with rationale.
- [x] 1.2 Remove `.circleci/config.yml` or replace it with an inactive/deprecated note.
- [x] 1.3 Move root `plan.md` to `docs/plans/archive/` or mark it historical at the top.

## Phase 2: Documentation Alignment

- [x] 2.1 Update `README.md` so support matrix, Quick Start, canonical API, and release package claims match current repo.
- [x] 2.2 Update `docs/ARCHITECTURE.md` to describe `+paraxial/` as canonical and `src/` as deprecated/transitional.
- [x] 2.3 Update `tests/README.md` to match `tests/test_all.m`, `portable_runner()`, and GitHub Actions behavior.
- [x] 2.4 Update `.atl/skill-registry.md` so Octave baseline and test command match README/CI.
- [x] 2.5 Create `docs/ROADMAP.md` with post-v2 cleanup phases: hygiene, docs, guardrails, packaging, legacy/addons.

## Phase 3: Guardrail Tests

- [x] 3.1 Create `tests/modern/test_RepositoryGuardrails.m` to scan `examples/canonical/*.m` for forbidden direct `src/beams` usage.
- [x] 3.2 Add README assertions: mentions `+paraxial/` as canonical and marks `src/` deprecated/transitional.
- [x] 3.3 Add BeamFactory registry assertions for supported beam names in `ParaxialBeams/BeamFactory.m`.
- [x] 3.4 Register `test_RepositoryGuardrails.m` in the existing portable test path (`tests/test_all.m` or `portable_runner()`).

## Phase 4: Packaging/Release Checklist

- [x] 4.1 Add release checklist to `docs/ROADMAP.md` or `docs/migration/` covering tests, version placeholder, artifacts, and changelog.
- [x] 4.2 Verify docs reference `.github/workflows/release.yml` rather than stale package instructions.
- [x] 4.3 Update `CHANGELOG.md` with an Unreleased cleanup section if changes are implemented.

## Phase 5: Verification

- [x] 5.1 Run `octave-cli --no-gui --eval "run('tests/test_all.m')"` and record result. Passed with 33 passed, 0 failed using 600s timeout.
- [x] 5.2 Run a docs/path audit: README, ARCHITECTURE, tests README, `.atl/skill-registry.md` agree on supported versions and commands.
- [x] 5.3 Confirm no beam physics files changed: `+paraxial/+beams/`, `src/beams/`, `src/parameters/`, `src/computation/` unchanged unless explicitly justified.
