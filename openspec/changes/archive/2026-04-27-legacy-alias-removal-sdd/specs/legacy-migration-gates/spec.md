# Spec: legacy-migration-gates

## Requirement 1 — Gate evidence traceability

The system SHALL maintain explicit evidence for Usage/Test/Docs/Release gates
before removing legacy aliases.

### Scenario: Verification report captures gate outcomes

- **GIVEN** a planned legacy alias removal change
- **WHEN** readiness verification is executed
- **THEN** a `verify-report.md` artifact records pass/fail for each gate
- **AND** includes references to command outputs or files used as evidence

## Requirement 2 — Guardrail enforcement

The system SHALL prevent reintroduction of deprecated aliases in canonical and
modern code surfaces.

### Scenario: Guardrail fails on deprecated alias usage

- **GIVEN** a file under `examples/canonical` or `tests/modern` references
  `HankeleHermite` or `HankeleLaguerre`
- **WHEN** `test_LegacyAliasGuardrail.m` runs
- **THEN** the test fails with a violation report
