# Archive Report: Performance Benchmarks SDD

## Change

- Change: `performance-benchmarks-sdd`
- Archived: 2026-04-27

## Specs Synced

| Domain | Action | Details |
|--------|--------|---------|
| `propagator-benchmark-harness` | No change | Delta only — no main spec created |
| `benchmark-baseline-output` | No change | Delta only |

Note: No main spec was promoted for this change.

## Archive Contents

- `proposal.md` ✅
- `design.md` ✅
- `tasks.md` ✅ — 12/12 tasks complete
- `verify-report.md` ✅ — verdict PASS
- `specs/propagator-benchmark-harness/spec.md` ✅
- `specs/benchmark-baseline-output/spec.md` ✅

## Verification Evidence

- Benchmark harness executed via GitHub Actions `workflow_dispatch` on PR #40.
- Baseline CSV captured: `docs/performance/baseline_quick_20260427_230128.csv`
- Quick mode results: gaussian_fft ~0.0108s, gaussian_analytic ~0.0056s, gaussian_raytrace ~14.03s (256x256 grid)
- Workflow run: `https://github.com/dreamjorge/Simulation_Scripts/actions/runs/25024076683`
- Artifact `benchmark-baseline` contains the CSV

## Source of Truth Updated

- `docs/performance/BENCHMARK_PROTOCOL.md` — protocol documentation
- `docs/performance/baseline_quick_20260427_230128.csv` — quick mode baseline
- `tests/performance/` — benchmark harness files
- `.github/workflows/octave.yml` — benchmark CI job added

## Notes

- Full benchmark mode not executed — only quick mode baseline captured. Quick mode sufficient per design.
- CI artifact publication (Phase 2 future work) not implemented — benchmark job does not fail CI on regression.