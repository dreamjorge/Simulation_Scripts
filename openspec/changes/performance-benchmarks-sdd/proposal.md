# Proposal: performance-benchmarks-sdd

## Why

After stabilizing architecture and legacy migration, we need quantitative
performance visibility across propagation strategies and representative problem
sizes.

Today, CI validates correctness but does not track runtime regressions. This
change introduces a lightweight, reproducible benchmark harness for Octave/
MATLAB workflows.

## What Changes

1. Add benchmark scripts for core propagators:
   - `FFTPropagator`
   - `AnalyticPropagator`
   - `RayTracePropagator`
2. Define controlled benchmark matrix (beam types, grid sizes, z-planes).
3. Emit machine-readable benchmark outputs (CSV/JSON) for trend tracking.
4. Add docs with execution protocol and interpretation guidance.

## Non-Goals

- No micro-optimizations in this change.
- No physics-model refactor.
- No hard CI fail thresholds yet (informational stage first).

## Success Criteria

- Benchmarks are reproducible with documented commands.
- Baseline runtime dataset is generated and committed (or stored artifact path).
- Developers can compare runs and detect regressions intentionally.
