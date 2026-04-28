# Package Distribution SDD — Tasks

## Phase 1: Local Package Scaffolding

### Task 1.1 — Create DESCR.ini
**File:** `DESCR.ini`
**Spec:** `octave-package-structure`

Create the Octave package metadata file at repo root.

Contents:
- `name=simulation_scripts`
- `version=` (placeholder, replaced by CI sed)
- `author=`, `maintainer=`, `description=` from project
- `homepage=` → GitHub repo URL
- `license=MIT` (or current LICENSE)
- `depends=` (empty — no required dependencies)

**Verification:** File is valid INI format; no Octave required.

---

### Task 1.2 — Create install.m
**File:** `install.m`
**Spec:** `octave-package-structure`, `matlab-toolbox-structure`

MATLAB/Octave install hook.

Behavior:
- Adds the repo root to the path (so `+paraxial/` namespaces are found)
- Adds all `+paraxial/**` subdirectories recursively
- Displays a message: `"Simulation_Scripts v<version> installed. See README for usage."`

Implementation approach:
- Use `mfilename('fullpath')` to find the install script's directory
- Use `addpath(genpath(...))` for MATLAB
- For Octave, the `pkg` tool calls `install.m` automatically after extracting

**Verification:** `install.m` runs without error in both Octave (`octave --eval "install"`) and MATLAB (`run install.m`).

---

### Task 1.3 — Create uninstall.m
**File:** `uninstall.m`
**Spec:** `octave-package-structure`, `matlab-toolbox-structure`

MATLAB/Octave uninstall hook.

Behavior:
- Removes the repo root from the path
- Removes all `+paraxial/**` subdirectories
- Displays a message confirming uninstall

**Verification:** `uninstall.m` runs without error; `which simulation_scripts_version` returns error after uninstall.

---

### Task 1.4 — Add simulation_scripts_version() to +paraxial/init.m
**File:** `+paraxial/init.m`
**Spec:** `version-reporting`

Modify or extend `+paraxial/init.m` to expose:

```matlab
function ver = simulation_scripts_version()
    [status, result] = system('git describe --tags --match "v*" --always');
    if status == 0
        ver = strtrim(result);
    else
        ver = '0.0.0-unknown';
    end
end
```

**Verification:**
- Octave: `octave --no-gui --eval "setpaths; ver = simulation_scripts_version(); disp(ver)"`
- MATLAB: After installing toolbox, `ver = simulation_scripts_version()` returns a string starting with `v`

---

## Phase 2: Release Workflow

### Task 2.1 — Create release.yml
**File:** `.github/workflows/release.yml`
**Spec:** `github-release-workflow`

Create the GitHub Actions workflow.

Trigger: `push:
  tags: ['v*']`

Jobs:
1. **test-and-build** (runs on `ubuntu-latest`)
   - Checkout with `fetch-depth: 0` (to get all tags)
   - Run `portable_runner.m` via `octave --no-gui --eval "setpaths; status = portable_runner(); if status ~= 0, error('tests failed'); end"`
   - Extract version from tag: `VERSION=${GITHUB_REF#refs/tags/v}`
   - Build Octave `.tar.gz` via `pkg tarball` (pkg is pre-installed in Ubuntu Octave)
   - Rename artifact: `mv simulation_scripts*.tar.gz simulation_scripts-${VERSION}.tar.gz`
   - Upload artifact: `actions/upload-release-artifact`
   - **MATLAB job** (runs on `ubuntu-latest` with MATLAB)
     - Same checkout and test step
     - Setup MATLAB via `matlab-actions/setup-matlab@v3` with `release: latest`
     - Build `.mltbx` via `matlab.addons.createToolbox(pwd, 'OutputFile', ['simulation_scripts-' version '-mltbx'])`
     - Upload artifact

2. **create-release** (runs on `ubuntu-latest`, needs `contents: write` permission)
   - Downloads both artifacts
   - Uses `softprops/action-gh-release@v2` or GitHub REST API to create a release
   - Attaches both artifacts
   - Auto-generates release notes via `generate_release_notes: true`

**Verification:** Tag `v2.0.0-test` locally, push, verify workflow triggers, check artifacts attached to draft release. Delete test tag after verification.

---

## Phase 3: README Update

### Task 3.1 — Update README.md installation section
**File:** `README.md`
**Spec:** `octave-package-structure`, `matlab-toolbox-structure`

Replace or supplement existing `setpaths.m` instructions with:

**Octave:**
```octave
pkg install 'https://github.com/dreamjorge/Simulation_Scripts/releases/latest/download/simulation_scripts-<VERSION>.tar.gz'
```

**MATLAB:**
```matlab
matlab.addons.install('https://github.com/dreamjorge/Simulation_Scripts/releases/latest/download/simulation_scripts-<VERSION>.mltbx')
```

Or double-click the `.mltbx` file.

Also add section on `simulation_scripts_version()` and link to the full CHANGELOG.

**Verification:** README renders correctly on GitHub; install commands are syntactically correct.

---

## Task Checklist

- [x] Task 1.1 — DESCR.ini created
- [x] Task 1.2 — install.m created and runs in Octave
- [x] Task 1.3 — uninstall.m created and runs in Octave
- [x] Task 1.4 — simulation_scripts_version() works in Octave
- [x] Task 2.1 — release.yml created and triggers on tag push
- [x] Task 3.1 — README.md installation section updated
