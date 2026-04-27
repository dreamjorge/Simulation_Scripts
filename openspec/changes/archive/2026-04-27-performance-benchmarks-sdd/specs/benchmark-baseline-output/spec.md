# Spec: benchmark-baseline-output

## Requirement 1 — Structured benchmark output

The benchmark run SHALL generate machine-readable output files.

### Scenario: CSV output generated

- **GIVEN** benchmark execution completes
- **WHEN** output is written
- **THEN** a CSV file includes scenario name, mode, grid size, repeats, and median runtime

## Requirement 2 — Baseline traceability

The project SHALL preserve a baseline dataset for future comparisons.

### Scenario: Baseline file recorded with run context

- **GIVEN** full benchmark mode is executed
- **WHEN** baseline is captured
- **THEN** baseline filename includes date/version context
- **AND** accompanying docs explain how to compare future runs to baseline
