# Archive Report: package-distribution-sdd

**Archived:** 2026-04-28
**Merge commit:** `fda7577` (PR #41 merged 2026-04-28T00:24:40Z)
**Origin:** `openspec/changes/package-distribution-sdd/`

## Summary

Created an installable package structure for Simulation_Scripts targeting both Octave and MATLAB, with automated GitHub release workflow.

## Scope Delivered

| Capability | Spec | Status |
|------------|------|--------|
| `octave-package-structure` | `specs/octave-package-structure/spec.md` | ✅ Implemented |
| `matlab-toolbox-structure` | `specs/matlab-toolbox-structure/spec.md` | ✅ Implemented |
| `version-reporting` | `specs/version-reporting/spec.md` | ✅ Implemented |
| `github-release-workflow` | `specs/github-release-workflow/spec.md` | ✅ Implemented |

## Files Committed

| File | Change |
|------|--------|
| `DESCR.ini` | New — Octave package metadata |
| `install.m` | New — Octave/MATLAB install hook |
| `uninstall.m` | New — Octave/MATLAB uninstall hook |
| `+paraxial/init.m` | Modified — added `simulation_scripts_version()` |
| `.github/workflows/release.yml` | New — release workflow |
| `README.md` | Modified — package installation instructions |

## Design Decisions

- **D1:** Git tags as version source of truth via `git describe --tags`
- **D2:** `simulation_scripts_version()` lives in `+paraxial/init.m` (not `src/version.m`)
- **D3:** Octave uses `pkg tarball` to create `.tar.gz`
- **D4:** MATLAB uses `matlab.addons.createToolbox` (no `package.xml` needed)
- **D5:** Release workflow runs `portable_runner.m` before packaging (fails fast on broken tests)

## Release Trigger

```bash
git tag v2.0.0 && git push origin v2.0.0
```

Workflow (`release.yml`) triggers on tag push, runs tests, builds `.tar.gz` (Octave) and `.mltbx` (MATLAB), and attaches both artifacts to the GitHub Release with auto-generated release notes.

## Verification

- Octave CI on PR #41: **passed** (1m50s)
- All tasks from `tasks.md` completed and checked off
- Manual install test: `install.m` runs without error in Octave
- `simulation_scripts_version()` returns git tag or `0.0.0-unknown`

## Active Changes After Archive

0 active SDD changes remain in `openspec/changes/`.