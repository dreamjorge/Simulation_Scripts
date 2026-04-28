# Package Distribution SDD — Design

## Context

Simulation_Scripts currently has no installable package structure. Users must manually call `setpaths.m` or add directories to their path. The goal is to make the project installable via `pkg install` (Octave) and `matlab.addons.install` (MATLAB), with versioning via Git tags and release automation via GitHub Actions.

## Design Decisions

### D1: Version source of truth — Git tags

The version string comes from the Git tag on the current commit, not a static file.

**Choice:** `git describe --tags --match "v*" --always` at runtime to compute the version string at install time.
**Rationale:** A `VERSION` file risks going out of sync with tags. Git is always authoritative. Both Octave and MATLAB can call `system("git describe ...")` or use `git2changelog` equivalents.
**Tradeoff:** Requires Git to be available at runtime — acceptable since both Octave and MATLAB run in Git-tracked environments (CI and typical user setups).

### D2: Version function lives in `+paraxial/init.m`

**Choice:** `simulation_scripts_version()` is implemented as a function in `+paraxial/init.m`.
**Rationale:** Octave and MATLAB both expose package version via the package root's `init.m` (Octave convention) or via the toolbox metadata (MATLAB convention). Having it in `+paraxial/init.m` means it is accessible immediately after the package namespace is on the path, with no additional files needed.
**Note:** `src/version.m` from the proposal is not created — `+paraxial/init.m` serves this purpose.

### D3: Octave package uses `pkg tarball` workflow

**Choice:** In CI, run `pkg tarball` on the repository root to produce `.tar.gz`.
**Rationale:** `pkg tarball` automatically includes all subdirectories. The `DESCR.ini` at the repo root provides metadata. The Octave package manager handles `install.m` automatically on `pkg install`.
**Layout:** `DESCR.ini` at repo root; `install.m` at repo root (called automatically by `pkg install`); `uninstall.m` at repo root (called by `pkg uninstall`).

### D4: MATLAB toolbox uses `matlab.addons.createToolbox`

**Choice:** In CI, use `matlab.addons.createToolbox(projectFolder, 'OutputFile', 'simulation_scripts-<VERSION>.mltbx')`.
**Rationale:** This is MATLAB's standard toolbox bundler. It bundles all folders in the project folder, including namespace packages (`+paraxial/`), into a single `.mltbx` file. `package.xml` is not needed — metadata is embedded in the toolbox properties.
**Note:** `package.xml` from the proposal is not created. The `package.xml` format is for MATLAB File Exchange, not for `.mltbx` bundles.

### D5: Release workflow runs tests before packaging

**Choice:** The release workflow runs `portable_runner.m` first. Only if it passes does it build and attach artifacts.
**Rationale:** Shipping broken packages is worse than not shipping. Attaching to an existing release on retag is acceptable for v2.0.0.

## File Layout After Implementation

```
Simulation_Scripts/          ← repo root
├── DESCR.ini                ← Octave package metadata
├── install.m                ← Octave/MATLAB install hook
├── uninstall.m              ← Octave/MATLAB uninstall hook
├── +paraxial/
│   └── init.m               ← package init + simulation_scripts_version()
├── .github/
│   └── workflows/
│       └── release.yml      ← release workflow
└── (existing files unchanged)
```

## Implementation Phases

### Phase 1 — Local package scaffolding
Files created locally (no CI needed to verify):
1. `DESCR.ini` — Octave package metadata
2. `install.m` — adds all `+paraxial/` namespaces to the path
3. `uninstall.m` — removes those paths
4. `+paraxial/init.m` — add `simulation_scripts_version()` function

### Phase 2 — CI / release workflow
`.github/workflows/release.yml` with:
1. Trigger on `refs/tags/v*`
2. Checkout tag
3. Run `portable_runner.m` via Octave (portable runner, no MATLAB needed)
4. Build `.tar.gz` via `pkg tarball`
5. Build `.mltbx` via `matlab.addons.createToolbox`
6. Attach both to GitHub Release
7. Auto-generate release notes

### Phase 3 — README update
Update `README.md` installation section to use the new package-based install instead of `setpaths.m`.

## Release Metadata (DESCR.ini)

```ini
[package]
name=simulation_scripts
version=__GET_FROM_GIT__
author=...
maintainer=...
description=Paraxial beam propagation and wavefront analysis in Octave and MATLAB
```

Version placeholder `__GET_FROM_GIT__` is replaced at CI build time via `sed`.

## CI Environment Variables

| Variable | Source | Purpose |
|---|---|---|
| `GITHUB_REF` | GitHub Actions | Tag name (e.g. `refs/tags/v2.0.0`) |
| `GITHUB_TOKEN` | GitHub Actions | Auth for release API |
| `VERSION` | `echo ${GITHUB_REF#refs/tags/v}` | Semantic version string |

## Key Risks and Mitigations

| Risk | Mitigation |
|---|---|
| `git describe` fails in CI artifact (detached HEAD) | Use `git fetch --tags` before describe; CI always has full tag context |
| MATLAB CI license issues | `matlab-actions/setup-matlab@v3` with `release: latest` uses GitHub's free MATLAB license for open source |
| `portable_runner.m` fails on tag | Workflow fails; no artifact published; release stays in draft |
| Octave `.tar.gz` exceeds GitHub release size limit (2GB) | Unlikely for this codebase; monitor first releases |

## Rollback

- If the release workflow fails mid-build: the GitHub Release remains in draft or is not published; retag after fix.
- If a published package is broken: users can download the previous tag's artifact directly from GitHub releases.
- `setpaths.m` remains functional as fallback for users who do not use the package install.
