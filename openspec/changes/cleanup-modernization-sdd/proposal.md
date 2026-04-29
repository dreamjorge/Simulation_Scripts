# Proposal: Cleanup & Modernization Alignment

## Intent

Resolve repository hygiene, documentation drift, and canonical API ambiguity after the v2 modernization work. The project already has `+paraxial/`, tests, CI, and release packaging; this change makes those decisions explicit, verifiable, and maintainable.

## Scope

### In Scope
- Remove or document stale local/tooling surfaces (`.circleci/`, `.opencode/`).
- Align public docs around supported Octave/MATLAB versions, canonical namespace, test commands, and release flow.
- Add architecture guardrails that verify canonical examples and public docs do not regress to deprecated surfaces.
- Harden packaging/release documentation without changing physical formulas.

### Out of Scope
- Rewriting beam physics, propagation algorithms, or parameter formulas.
- Removing `src/` or legacy examples.
- Changing package versioning semantics beyond documentation/checklists.
- Running MATLAB/Octave builds as part of this planning change.

## Capabilities

### New Capabilities
- `repository-modernization-hygiene`: Rules for repository hygiene, stale CI/tooling removal, and documentation consistency.
- `canonical-api-guardrails`: Structural checks that keep canonical API usage on `+paraxial/` and `BeamFactory`.

### Modified Capabilities
- `ci-test-status`: CI/test documentation SHALL match the actual GitHub Actions and runner behavior.
- `paraxial-beam-factory`: Public docs SHALL describe `BeamFactory` and `+paraxial/` as the canonical API surface.

## Approach

Use a documentation-first cleanup with small guardrail tests. First normalize repository metadata and docs, then add low-cost structural tests to prevent drift. Keep implementation away from numerical behavior.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `.circleci/config.yml` | Removed/Archived | Stale welcome-only CI no longer matches GitHub Actions reality. |
| `.opencode/` / `.gitignore` | Modified | Decide and encode local tooling policy. |
| `README.md`, `docs/ARCHITECTURE.md`, `tests/README.md`, `.atl/skill-registry.md` | Modified | Align support matrix, commands, and canonical API story. |
| `tests/modern/` | Modified | Add structural guardrail tests. |
| `docs/ROADMAP.md` | Created | Current post-v2 modernization roadmap. |
| `.github/workflows/*.yml` | Reference | Ensure docs reflect active CI jobs. |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Removing stale config surprises users | Low | Document GitHub Actions as canonical before deleting CircleCI. |
| Guardrail tests become brittle | Med | Test stable invariants only: paths, docs, canonical examples. |
| Docs overpromise package behavior | Med | Verify claims against existing files before editing. |

## Rollback Plan

Revert the cleanup commit(s). Restore `.circleci/config.yml` from git history if required. Remove new guardrail tests if they block valid future migrations.

## Dependencies

- Existing GitHub Actions workflows: `.github/workflows/octave.yml`, `matlab.yml`, `release.yml`.
- Existing SDD specs: `ci-test-status`, `paraxial-beam-factory`.

## Success Criteria

- [ ] No stale CI/tooling config is left undocumented.
- [ ] Docs agree on supported versions, test commands, canonical API, and release flow.
- [ ] Canonical examples avoid deprecated `src/` direct usage.
- [ ] New guardrail tests are registered in the portable suite.
- [ ] No beam formulas or numerical APIs are changed.
