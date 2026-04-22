# Spec: legacy-alias-removal

## Requirement 1 — Controlled alias removal

The system SHALL remove deprecated `Hankele*` aliases only after readiness
gates are verified and recorded.

### Scenario: Alias files removed after gate pass

- **GIVEN** all readiness gates are marked passed in verification artifacts
- **WHEN** the alias removal change is applied
- **THEN** `legacy/compat/HankeleHermite.m` and `legacy/compat/HankeleLaguerre.m`
  are removed from the repository
- **AND** a release checkpoint records the commit hash and rollback instructions

## Requirement 2 — Modern API stability

The system SHALL preserve canonical and modern API behavior during alias removal.

### Scenario: Portable suite remains green

- **GIVEN** alias removal patch is applied
- **WHEN** `tests/portable_runner.m` executes
- **THEN** test run finishes with zero failures
