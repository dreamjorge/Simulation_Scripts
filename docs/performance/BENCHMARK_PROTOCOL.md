# Benchmark Protocol (Propagators)

This document defines how to run and interpret the benchmark harness added in `tests/performance/`.

## Scope

Benchmarks are **informational** (not pass/fail correctness checks). They measure runtime trends for:

- `FFTPropagator`
- `AnalyticPropagator`
- `RayTracePropagator`

Scenarios include Gaussian workloads plus representative Hermite/Laguerre medium-grid cases.

## Files

- `tests/performance/bench_config.m` — benchmark matrix and mode definitions (`quick`/`full`)
- `tests/performance/bench_utils.m` — timing/CSV helper utilities
- `tests/performance/run_benchmarks.m` — benchmark runner entrypoint

## Modes

- **quick**: low-latency local sanity pass
  - warm-up: 1
  - repeats: 3
  - includes FFT + Analytic + RayTrace Gaussian scenarios

- **full**: baseline capture mode
  - warm-up: 2
  - repeats: 5
  - full matrix across grid tiers + representative beam families

## Output

CSV files are written to:

`docs/performance/`

Schema:

- `timestamp`
- `mode`
- `scenario`
- `propagator`
- `beam_type`
- `grid_tier`
- `nx`, `ny`
- `z_final`
- `warmup`
- `repeats`
- `median_runtime_sec`
- `min_runtime_sec`
- `max_runtime_sec`

When running `full` mode, an extra baseline file is written using date context:

- `baseline_YYYY-MM-DD.csv`

## Commands

### GNU Octave

```bash
octave --quiet --eval "cd('tests/performance'); run_benchmarks('quick');"
octave --quiet --eval "cd('tests/performance'); run_benchmarks('full');"
```

### MATLAB

```bash
matlab -batch "cd('tests/performance'); run_benchmarks('quick');"
matlab -batch "cd('tests/performance'); run_benchmarks('full');"
```

## Interpretation Rules

1. Compare **median** runtimes, not single samples.
2. Compare runs generated on similar machine/software environments.
3. Treat large regressions as investigation triggers, not immediate failures.
4. Prefer trend analysis per scenario (`small -> medium -> large`) over absolute numbers.

## CI Follow-up (future phase)

Proposed next step:

- publish benchmark CSV as CI artifacts
- keep benchmark stage non-blocking
- optionally warn on configurable regression thresholds
