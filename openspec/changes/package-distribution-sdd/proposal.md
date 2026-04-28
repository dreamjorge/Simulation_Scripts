# Proposal: Package Distribution SDD

## Intent

Formalize distribution of Simulation_Scripts as an installable package for both Octave and MATLAB, enabling versioning, registry-free installation, and GitHub-based release management.

## Scope

### In Scope
- Create Octave package structure (`DESCR.ini`, `.tar.gz` packaging)
- Create MATLAB toolbox structure (`.mltbx` bundle)
- Add version management (semantic versioning via Git tags, `VERSION` file)
- Create `install.m` / `uninstall.m` scripts for MATLAB
- Create `install.m` / `uninstall.m` scripts for Octave
- Update `README.md` with installation instructions
- Add GitHub Actions release workflow to build and attach packages to releases

### Out of Scope
- Publishing to MATLAB File Exchange or Octave Forge (requires account approval)
- Dependency resolution between packages
- Supporting MATLAB < R2020b or Octave < 11.1.0

## Capabilities

### New Capabilities
- `octave-package-structure`: Octave package format with DESCR.ini, installation scripts, and .tar.gz distribution
- `matlab-toolbox-structure`: MATLAB toolbox (.mltbx) format with installation scripts
- `version-reporting`: `ver = simulation_scripts_version()` returning semantic version string
- `github-release-workflow`: CI workflow to build packages and attach to GitHub releases on tag

### Modified Capabilities
- None — this is net-new infrastructure, not a change to existing functionality

## Approach

**Octave**: Use the standard Octave package format. Package root contains `DESCR.ini` (metadata), `index.ttk` (if needed), and the namespace folders `+paraxial/`, `+paraxial/+beams/`, etc. Create `.tar.gz` archive via `pkg` tool in CI. Installation: `pkg install simulation_scripts-<version>.tar.gz`.

**MATLAB**: Use Toolbox format (.mltbx). Bundle all `+paraxial/` namespace folders, `ParaxialBeams/`, `src/`, and `tests/` into a single toolbox file. Installation: double-click `.mltbx` or `matlab.addons.install()`.

**Versioning**: Use Git tags (`v2.0.0`, `v2.1.0`, etc.) to drive semantic version. Create `src/version.m` or `+paraxial/init.m` to expose `ver = simulation_scripts_version()`.

**GitHub Actions**: On tag push (`refs/tags/v*`), run a job that:
1. Builds `.tar.gz` (Octave) and `.mltbx` (MATLAB) from current source
2. Attaches artifacts to the GitHub Release

**Installation scripts**: `install.m` adds paths; `uninstall.m` removes them. For Octave, `install.m` is called by the package manager automatically.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `+paraxial/init.m` | Modified | Add version reporting function |
| `ParaxialBeams/setpaths.m` | Modified | Keep as fallback; document install-based approach |
| `src/version.m` | New | Version reporting |
| `install.m` (root) | New | MATLAB/Octave install script |
| `uninstall.m` (root) | New | MATLAB/Octave uninstall script |
| `DESCR.ini` (root) | New | Octave package metadata |
| `package.xml` (root) | New | MATLAB toolbox metadata |
| `.github/workflows/release.yml` | New | Package build and release workflow |
| `README.md` | Modified | Update installation section |
| `CHANGELOG.md` | Modified | Ensure format supports release versioning |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Octave Forge vs GitHub releases strategy unclear | Medium | Use GitHub releases as primary; note Forge as optional future |
| MATLAB toolbox build requires MATLAB license in CI | High | Use `matlab-actions/setup-matlab@v3` with `release: latest` — free for open source |
| Package contents differ between Octave/MATLAB | Low | Both bundle same source tree; CI validates with portable_runner |

## Rollback Plan

- If release workflow fails: do not tag; revert tag push
- If package is broken: previous release remains available on GitHub; users downgrade via `pkg install` URL or manual download
- Installation scripts are additive — existing `setpaths.m` remains functional fallback

## Dependencies

- GitHub Actions with `matlab-actions/setup-matlab@v3` (free for open source repos via GitHub's MATLAB action)
- Octave `pkg` tool for .tar.gz creation (available in standard Ubuntu Octave CI image)
- Semantic versioning discipline (tags must follow `vM.m.p` format)

## Success Criteria

- [ ] `octave --no-gui --eval "pkg install simulation_scripts-*.tar.gz"` installs all namespaces
- [ ] `matlab.addons.install('simulation_scripts.mltbx')` installs toolbox without manual addpath
- [ ] `ver = simulation_scripts_version()` returns `'2.0.0'` (or current version)
- [ ] GitHub release on tag `v2.0.0` contains both `.tar.gz` and `.mltbx` artifacts
- [ ] `tests/portable_runner.m` passes after package install (no manual path setup)
- [ ] Uninstall cleanly removes all added paths