# Spec: propagator-benchmark-harness

## Requirement 1 — Reproducible benchmark runner

The system SHALL provide a deterministic benchmark runner for core propagators.

### Scenario: Run benchmark matrix in quick mode

- **GIVEN** the benchmark harness is present
- **WHEN** a developer runs quick mode
- **THEN** FFT, Analytic, and RayTrace benchmark scenarios execute
- **AND** the run completes with structured results for each scenario

## Requirement 2 — Stable timing protocol

The system SHALL use warm-up and repeated measurements before reporting.

### Scenario: Median runtime is reported per scenario

- **GIVEN** a scenario in the benchmark matrix
- **WHEN** runner executes timing loop
- **THEN** warm-up iterations are excluded from final metric
- **AND** reported metric uses median of repeated measured runs
