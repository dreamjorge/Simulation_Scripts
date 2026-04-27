# Design: performance-benchmarks-sdd

## Context

Correctness tests are already comprehensive. We now need a dedicated benchmark
layer to answer:

- Which propagator is fastest for each workload shape?
- How runtime scales with grid size and z-plane count?
- Whether recent refactors changed runtime behavior materially?

## Approach

Use a deterministic benchmark runner under `tests/performance/` with:

- warm-up phase (avoid first-run noise)
- fixed random seeds where applicable
- repeated runs + median runtime
- output as structured files (`.csv` and optional `.json`)

### Candidate benchmark matrix

1. Gaussian + FFT propagation
2. Gaussian + Analytic propagation
3. Gaussian + RayTrace propagation
4. Hermite/Laguerre representative cases (medium grid)

Grid tiers:

- small: 256x256
- medium: 512x512
- large: 1024x1024

## Deliverables

- `tests/performance/run_benchmarks.m`
- `tests/performance/bench_config.m`
- `tests/performance/bench_utils.m`
- `docs/performance/BENCHMARK_PROTOCOL.md`
- baseline output file (e.g., `docs/performance/baseline_2026-04-22.csv`)

## Risks and Mitigations

- **Risk:** noisy timings across machines.
  - **Mitigation:** document machine/context; compare trends, not absolutes.

- **Risk:** benchmark runtime too long for local dev.
  - **Mitigation:** add quick mode (small matrix) and full mode.

## Rollout

Phase 1: informational benchmark tooling + baseline generation.

Phase 2 (future): optional CI artifact publication and soft regression alerts.
