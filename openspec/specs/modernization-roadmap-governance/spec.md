# Modernization Roadmap Governance Specification

## Purpose

Keep active modernization work discoverable, scoped, and aligned with executable guardrails.

## Requirements

### Requirement: Roadmap Owns Active Next Steps

`docs/ROADMAP.md` MUST identify active cleanup phases and SHOULD point to current SDD changes for detailed execution.

#### Scenario: New modernization SDD exists

- GIVEN an active OpenSpec change exists for modernization cleanup
- WHEN `docs/ROADMAP.md` is reviewed
- THEN it SHOULD reference the change or its tracked work area
- AND root `plan.md` MUST NOT be presented as the active roadmap.

### Requirement: Compatibility Reduction Is Planned Separately

Reducing or removing deprecated `src/` behavior MUST require a dedicated plan with preconditions.

#### Scenario: `src/` cleanup is discussed

- GIVEN docs mention future `src/` reduction
- WHEN the work is scoped
- THEN docs MUST state that `src/` removal is out of scope for cleanup-only changes
- AND preconditions MUST include compatibility tests and user-facing migration notes.

### Requirement: Guardrails Track Stable Invariants

Guardrail tests SHOULD verify stable architecture invariants, not exact prose formatting.

#### Scenario: Docs wording changes

- GIVEN docs are edited without changing policy
- WHEN guardrail tests run
- THEN tests SHOULD pass if canonical CI, canonical API, and deprecated-surface policy remain present
- AND tests SHOULD fail if those policies disappear.
