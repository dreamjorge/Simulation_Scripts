# Proposal: Post-v2 Modernization Next Steps

## Intent

Turn the remaining cleanup roadmap into executable SDD artifacts. The repo already completed a first hygiene pass; this change defines the next bounded wave: runner path clarity, docs/guardrail alignment, addons inventory, and a future compatibility reduction plan without touching beam physics.

## Scope

### In Scope
- Clarify canonical vs deprecated path setup in `tests/portable_runner.m` and related docs.
- Add/adjust guardrails for runner path policy and roadmap/doc consistency.
- Inventory `ParaxialBeams/Addons/` into runtime, plotting, vendored, or removable buckets.
- Draft a compatibility reduction plan for `src/` without removing `src/`.

### Out of Scope
- Rewriting optical formulas, propagation algorithms, or numerical behavior.
- Deleting `src/`, legacy examples, or addons without a dedicated follow-up change.
- Building release/package artifacts.
- Changing public `BeamFactory.supportedTypes()` names.

## Capabilities

### New Capabilities
- `runner-path-policy`: Canonical and deprecated MATLAB/Octave paths are explicit and testable.
- `legacy-addons-inventory`: Addons are classified before removal or migration decisions.
- `modernization-roadmap-governance`: Active roadmap, SDD changes, and docs remain aligned.

### Modified Capabilities
- None.

## Approach

Use a documentation-and-guardrail-first cleanup. Update runner setup comments/structure, encode stable invariants in guardrail tests, inventory addons without deletion, and document future migration gates.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `tests/portable_runner.m` | Modified | Separate canonical, deprecated, utility, and legacy paths. |
| `tests/modern/test_RepositoryGuardrails.m` | Modified | Add runner/path and roadmap consistency assertions. |
| `docs/ROADMAP.md` | Modified | Track next-wave cleanup and compatibility reduction gates. |
| `docs/ADDONS_INVENTORY.md` | Created | Classify addon files with rationale and follow-up action. |
| `docs/COMPATIBILITY_REDUCTION.md` | Created | Define preconditions for reducing deprecated `src/` usage. |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Path cleanup breaks legacy tests | Medium | Preserve `src/` paths but label them deprecated compatibility. |
| Addons classification is incomplete | Medium | Inventory only; no deletions in this change. |
| Guardrails become brittle | Medium | Assert stable policy text and paths, not prose formatting. |

## Rollback Plan

Revert the docs/guardrail/runner commits. Since no runtime formulas or deletions are included, rollback restores previous path setup and documentation only.

## Dependencies

- Existing roadmap: `docs/ROADMAP.md`.
- Existing guardrails: `tests/modern/test_RepositoryGuardrails.m`.
- Canonical test runner: `tests/portable_runner.m`.

## Success Criteria

- [ ] Runner setup explicitly adds repo root for `+paraxial/` and labels `src/` as deprecated compatibility.
- [ ] Guardrails verify runner/docs do not regress.
- [ ] Addons inventory exists with every top-level addon classified.
- [ ] Compatibility reduction plan exists and does not remove `src/`.
- [ ] Portable Octave runner passes after implementation.
