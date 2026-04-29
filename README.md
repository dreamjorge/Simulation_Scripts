# Simulation_Scripts

Scripts for simulation of optical beam propagation (Gaussian, Hermite-Gauss, Laguerre-Gauss beams).

**Author:** Ugalde-Ontiveros J.A.

[![Octave CI](https://github.com/dreamjorge/Simulation_Scripts/actions/workflows/octave.yml/badge.svg)](https://github.com/dreamjorge/Simulation_Scripts/actions/workflows/octave.yml)
[![MATLAB CI](https://github.com/dreamjorge/Simulation_Scripts/actions/workflows/matlab.yml/badge.svg)](https://github.com/dreamjorge/Simulation_Scripts/actions/workflows/matlab.yml)

## Project Status

GitHub Actions is the canonical CI system for this repository. The active workflows are:

- `.github/workflows/octave.yml` — Octave portable, legacy compatibility, and benchmark checks.
- `.github/workflows/matlab.yml` — MATLAB portable and legacy compatibility checks.
- `.github/workflows/release.yml` — release artifact generation for tagged versions.

The current modernization roadmap lives in `docs/ROADMAP.md`; historical implementation plans remain under `docs/plans/` or are marked as historical.

## Project Structure

```
Simulation_Scripts/
├── +paraxial/                  % Canonical package namespace for new code
│   ├── +beams/                 % Beam classes
│   ├── +parameters/            % Beam parameter classes
│   ├── +computation/           % Formula/logic layer
│   ├── +propagation/           % Field and ray propagation strategies
│   └── +visualization/         % Visualization and wavefront utilities
├── src/                        % Deprecated transition adapters / legacy layout
│   ├── beams/                  % Beam classes
│   │   ├── ParaxialBeam.m      % ⭐ Abstract base class
│   │   ├── GaussianBeam.m
│   │   ├── HermiteBeam.m
│   │   ├── LaguerreBeam.m
│   │   ├── ElegantHermiteBeam.m
│   │   ├── ElegantLaguerreBeam.m
│   │   ├── HankelHermite.m
│   │   ├── HankelLaguerre.m
│   │   └── (modern classes only)
│   ├── parameters/             % Beam parameter classes
│   │   ├── GaussianParameters.m
│   │   ├── HermiteParameters.m
│   │   ├── LaguerreParameters.m
│   │   ├── ElegantHermiteParameters.m
│   │   └── ElegantLaguerreParameters.m
│   ├── propagation/
│   │   ├── field/              % Field-based propagation
│   │   │   ├── IPropagator.m   % ⭐ Strategy interface
│   │   │   ├── FFTPropagator.m
│   │   │   └── AnalyticPropagator.m
│   │   └── rays/                % Ray-based propagation
│   │       ├── RayTracePropagator.m
│   │       ├── RayBundle.m
│   │       ├── RayTracer.m
│   │       ├── OpticalRay.m
│   │       └── CylindricalRay.m
│   └── visualization/
│       └── VisualizationUtils.m
├── ParaxialBeams/              % Utilities
│   ├── PhysicalConstants.m
│   ├── GridUtils.m
│   ├── FFTUtils.m
│   ├── AnalysisUtils.m
│   ├── PolynomialUtils.m
│   ├── BeamFactory.m
│   └── Addons/
├── examples/
│   ├── canonical/          % ✅ Canonical examples for new users
│   │   ├── MainGauss_refactored.m
│   │   ├── MainMultiMode.m
│   │   ├── ExampleRayTracing.m
│   │   └── ExampleHankelPropagation.m
│   └── legacy/
│       ├── archive/        % Old API examples (for reference)
│       ├── generators/     % Figure generators for papers
│       ├── research/       % Thesis-specific scripts
│       └── LEGACY_POLICY.md
├── tests/                  % Portable Octave/MATLAB test suite
├── setpaths.m              % Path initialization utility
├── legacy/
│   └── compat/             % Compatibility docs (Hankele* aliases removed)
├── docs/
│   ├── ARCHITECTURE.md    % Architecture documentation
│   └── ROADMAP.md         % Active modernization roadmap
└── README.md
```

## Beam API Contract

Every beam type inherits from `ParaxialBeam` and implements:

| Method | Returns | Description |
|--------|---------|-------------|
| `opticalField(X, Y, z)` | `[Ny x Nx] complex` | Complex field on Cartesian grid |
| `getParameters(z)` | `GaussianParameters` | Beam params at position z |
| `beamName()` | `char` | String identifier like 'hermite_3_2' |

## Namespace Convention

**`+paraxial/` is the canonical namespace** for beam classes. The `src/` directory is deprecated but remains functional during the Strangler Fig migration transition.

| Location | Status | Usage |
|----------|--------|-------|
| `+paraxial/+beams/*.m` | ✅ Canonical | New code should use these classes directly |
| `src/beams/*.m` | ⚠️ Deprecated | Emits warning on instantiation; functional but being phased out |
| `BeamFactory.create()` | ✅ Preferred | Routes to `+paraxial/` automatically |

## Quick Start

### Install the Package

**Octave:**

```matlab
pkg install 'https://github.com/dreamjorge/Simulation_Scripts/releases/latest/download/simulation_scripts-<VERSION>.tar.gz'
```

Or download a specific release from the [GitHub Releases](https://github.com/dreamjorge/Simulation_Scripts/releases) page and install locally:

```matlab
pkg install simulation_scripts-<VERSION>.tar.gz
```

**MATLAB:**

Double-click the `.mltbx` file from the [GitHub Releases](https://github.com/dreamjorge/Simulation_Scripts/releases) page, or install from URL:

```matlab
matlab.addons.install('https://github.com/dreamjorge/Simulation_Scripts/releases/latest/download/simulation_scripts-<VERSION>.mltbx')
```

> **Note:** Replace `<VERSION>` with the desired release version (e.g. `2.0.0`). Omit the `v` prefix.

To check the installed version:

```matlab
ver = simulation_scripts_version()
```

### Manual Installation (fallback)

If you prefer not to use the package, add paths manually:

```matlab
% Option 1: Use setpaths() utility (adds +paraxial/, utilities, and adapters)
setpaths

% Option 2: Use +paraxial/ directly (recommended for direct class access)
addpath('+paraxial/+beams', '+paraxial/+parameters', '+paraxial/+computation');
addpath('ParaxialBeams');

% Create beam via Factory (routes to +paraxial/ automatically)
beam = BeamFactory.create('gaussian', 100e-6, 632.8e-9);

% Use beam directly
grid = GridUtils(1024, 1024, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
field = beam.opticalField(X, Y, 0);

% Or propagate via FFT
prop = FFTPropagator(grid, 632.8e-9);
field_at_z = prop.propagate(beam, 0.1);
```

## Design Patterns

### Strategy Pattern — Propagators

Three interchangeable propagation methods:

```matlab
% FFT (angular spectrum)
prop = FFTPropagator(grid, lambda);
field = prop.propagate(beam, z);

% Analytic (direct formula)
prop = AnalyticPropagator(grid);
field = prop.propagate(beam, z);

% Ray tracing
prop = RayTracePropagator(grid, 'RK4', 1e-3);
bundle = prop.propagate(beam, z);
```

### Factory Pattern — BeamFactory

```matlab
% All beams via Factory
g  = BeamFactory.create('gaussian', 100e-6, 632.8e-9);
hg = BeamFactory.create('hermite', 100e-6, 632.8e-9, 'n', 2, 'm', 1);
lg = BeamFactory.create('laguerre', 100e-6, 632.8e-9, 'l', 1, 'p', 0);
```

## Canonical Examples

Recommended examples for new users (in `examples/canonical/`):

| File | Description |
|------|-------------|
| `examples/canonical/MainGauss_refactored.m` | Gaussian beam propagation |
| `examples/canonical/MainMultiMode.m` | Multi-mode Hermite/Laguerre |
| `examples/canonical/ExampleRayTracing.m` | Ray tracing visualization |

## Constants and Utilities

### PhysicalConstants

```matlab
PC = PhysicalConstants;
k = PC.waveNumber(lambda);
zr = PC.rayleighDistance(w0, lambda);
R = PC.radiusOfCurvature(z, zr);
gouy = PC.gouyPhase(z, zr);
```

### GridUtils

```matlab
grid = GridUtils(Nx, Ny, Dx, Dy);
[X, Y] = grid.create2DGrid();
[Kx, Ky] = grid.createFreqGrid();
[r, theta] = grid.createPolarGrid();
```

### FFTUtils

```matlab
fftOps = FFTUtils(true, true);  % normalize, shift
G = fftOps.fft2(g);
g = fftOps.ifft2(G);
```

## Compatibility

- **GNU Octave 11.1.0+**
- **MATLAB R2020b+**

No `classdef` folders are used. All files are individual `.m` files.

## Tests

`tests/test_all.m` is a thin wrapper around `tests/portable_runner.m`, which is the canonical non-interactive runner used by CI.

```bash
# Octave
octave --no-gui --eval "run('tests/test_all.m')"

# Octave (legacy-only fast check)
octave --no-gui --eval "run('tests/legacy_compat/run_legacy_compat.m')"

# MATLAB
matlab -batch "run('tests/test_all.m')"
```

CI uses the same runner and fails the job when `portable_runner()` returns a non-zero failure count.

## Version

The package version is derived from Git tags. To check the version programmatically:

```matlab
ver = simulation_scripts_version()
```

This returns the Git tag (e.g. `'v2.0.0'`) or `'0.0.0-unknown'` if Git is not available.

## Release Packaging

Tagged releases (`v*`) trigger `.github/workflows/release.yml`, which runs the portable tests and uploads:

- `simulation_scripts-<VERSION>.tar.gz` for GNU Octave packages.
- `simulation_scripts-<VERSION>.mltbx` for MATLAB toolbox installs.

Before tagging a release, follow the release checklist in `docs/ROADMAP.md`.

## References

- Kogelnik, H., & Li, T. (1966). Laser beams and resonators. *Applied Optics*.
- Siegman, A. E. (1986). *Lasers*. University Science Books.
- Detailed references by implementation (Gaussian/Hermite/Laguerre/Elegant/Hankel): `docs/ARCHITECTURE.md` -> `References by Implementation`.

## Legacy Migration

- Incremental plan (Strangler pattern): `docs/migration/LEGACY_MIGRATION_PLAN.md`
- Release checkpoint (2026-04-15): `docs/migration/RELEASE_CHECKPOINT_2026-04-15.md`
- Release checkpoint (2026-04-22): `docs/migration/RELEASE_CHECKPOINT_2026-04-22.md`
- Alias removal release plan: `docs/migration/ALIAS_REMOVAL_RELEASE_PLAN.md`

## Uninstall

**Octave:**

```matlab
pkg uninstall simulation_scripts
```

**MATLAB:**

```matlab
matlab.addons.uninstall('Simulation_Scripts')
```

Or run `uninstall.m` manually.

## Changelog

- `CHANGELOG.md`
