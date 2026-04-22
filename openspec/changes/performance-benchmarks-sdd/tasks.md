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

- [ ] 3.1 Execute quick benchmark mode and validate output schema
- [ ] 3.2 Execute full benchmark mode and capture baseline dataset
- [ ] 3.3 Sanity-check benchmark trends vs expected complexity

> Pending due to local runtime constraints in this environment:
> - `octave` executable is unavailable in PATH
> - MATLAB CLI fails license checkout (License Manager Error -96)

## Phase 4: Documentation and handoff

- [x] 4.1 Add `docs/performance/BENCHMARK_PROTOCOL.md`
- [x] 4.2 Document command snippets for Octave/MATLAB
- [x] 4.3 Define follow-up for CI artifact integration (future phase)
