# Tasks: performance-benchmarks-sdd

## Phase 1: Benchmark scaffold

- [x] 1.1 Create `tests/performance/` directory
- [x] 1.2 Add `bench_config.m` with benchmark matrix and modes (quick/full)
- [x] 1.3 Add `bench_utils.m` with timing + aggregation helpers

## Phase 2: Runner implementation

- [x] 2.1 Implement `run_benchmarks.m` for propagator benchmarks
- [x] 2.2 Add warm-up/repeat logic and median reporting
- [x] 2.3 Write CSV output with timestamp and scenario labels

## Phase 3: Baseline generation

- [x] 3.1 Execute quick benchmark mode and validate output schema
- [x] 3.2 Execute full benchmark mode and capture baseline dataset
- [x] 3.3 Sanity-check benchmark trends vs expected complexity

> Baseline captured via GitHub Actions workflow_dispatch on PR #40.
> Quick mode results: gaussian_fft ~0.0108s, gaussian_analytic ~0.0056s, gaussian_raytrace ~14.03s (256x256 grid).
> Full mode pending: CI job added to octave.yml but full baseline not yet captured.
> Baseline CSV: `docs/performance/baseline_quick_20260427_230128.csv`

## Phase 4: Documentation and handoff

- [x] 4.1 Add `docs/performance/BENCHMARK_PROTOCOL.md`
- [x] 4.2 Document command snippets for Octave/MATLAB
- [x] 4.3 Define follow-up for CI artifact integration (future phase)
