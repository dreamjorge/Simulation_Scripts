# Tasks: Normalize Wavefront Test Paths

## Phase 1: Path Setup Correction

- [x] 1.1 Update `tests/modern/test_Wavefront.m` to derive repo root from `tests/modern/../..`.
- [x] 1.2 Add repo root to the path for `+paraxial/` package parent semantics.
- [x] 1.3 Keep `src/*`, `ParaxialBeams/`, and `ParaxialBeams/Addons/` additions rooted at the actual repo root.

## Phase 2: Guardrails

- [x] 2.1 Extend `tests/modern/test_RepositoryGuardrails.m` to detect stale `repoRoot = fullfile(testDir, '..')` modern-test path setup when project paths are added from it.
- [x] 2.2 Ensure the guardrail remains focused on stable path semantics, not exact formatting.

## Phase 3: Verification

- [x] 3.1 Run direct Wavefront test with Octave CLI and confirm no missing-path `addpath` warnings. Passed; only expected deprecated `src/beams/GaussianBeam` warning remains.
- [x] 3.2 Run direct repository guardrail test and confirm all guardrails pass. Passed with 14/14 guardrails.
- [x] 3.3 Run full portable Octave suite and confirm 0 failures. Passed with 33 passed, 0 failed.

## Phase 4: OpenSpec Closure Prep

- [x] 4.1 Record verification results in the eventual apply/verify report.
- [x] 4.2 Confirm no beam physics files changed under `+paraxial/+beams/`, `src/beams/`, `src/parameters/`, or `src/computation/`.
