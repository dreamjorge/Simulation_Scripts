# Simulation_Scripts

Scripts for simulation of optical beam propagation (Gaussian, Hermite-Gauss, Laguerre-Gauss beams).

**Author:** Ugalde-Ontiveros J.A.

## Estructura

```
Simulation_Scripts/
├── src/                        % Modern library (organized by responsibility)
│   ├── beams/                  % Beam classes
│   │   ├── ParaxialBeam.m      % ⭐ Abstract base class
│   │   ├── GaussianBeam.m
│   │   ├── HermiteBeam.m
│   │   ├── LaguerreBeam.m
│   │   ├── ElegantHermiteBeam.m
│   │   ├── ElegantLaguerreBeam.m
│   │   ├── HankelHermite.m
│   │   ├── HankelLaguerre.m
│   │   └── HankeleHermite.m, HankeleLaguerre.m (legacy aliases)
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
│   │   └── ExampleRayTracing.m
│   └── ... (legacy examples)
├── tests/                  % Test suite (~380 tests)
├── setpaths.m              % Path initialization utility
├── docs/
│   └── ARCHITECTURE.md    % Architecture documentation
└── README.md
```

## Beam API Contract

Every beam type inherits from `ParaxialBeam` and implements:

| Method | Returns | Description |
|--------|---------|-------------|
| `opticalField(X, Y, z)` | `[Ny x Nx] complex` | Complex field on Cartesian grid |
| `getParameters(z)` | `GaussianParameters` | Beam params at position z |
| `beamName()` | `char` | String identifier like 'hermite_3_2' |

## Uso Rápido

### MATLAB/Octave

```matlab
% Option 1: Use setpaths() utility
setpaths

% Option 2: Add paths manually
addpath('src/beams', 'src/parameters', 'src/propagation/field', 'src/propagation/rays', 'src/visualization');
addpath('ParaxialBeams', 'ParaxialBeams/Addons');

% Crear beam via Factory
beam = BeamFactory.create('gaussian', 100e-6, 632.8e-9);

% Usar beam directamente
grid = GridUtils(1024, 1024, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
field = beam.opticalField(X, Y, 0);

% O propagar via FFT
prop = FFTPropagator(grid, 632.8e-9);
field_at_z = prop.propagate(beam, 0.1);
```

## Patrones de Diseño

### Strategy Pattern — Propagators

Tres métodos de propagación intercambiables:

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
% Todos los beams vía Factory
g  = BeamFactory.create('gaussian', 100e-6, 632.8e-9);
hg = BeamFactory.create('hermite', 100e-6, 632.8e-9, 'n', 2, 'm', 1);
lg = BeamFactory.create('laguerre', 100e-6, 632.8e-9, 'l', 1, 'p', 0);
```

## Canonical Examples

Ejemplos recomendados para nuevos usuarios (en `examples/canonical/`):

| File | Description |
|------|-------------|
| `examples/canonical/MainGauss_refactored.m` | Gaussian beam propagation |
| `examples/canonical/MainMultiMode.m` | Multi-mode Hermite/Laguerre |
| `examples/canonical/ExampleRayTracing.m` | Ray tracing visualization |

## Constantes y Utilidades

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

## Compatibilidad

- **GNU Octave 11.1.0+**
- **MATLAB R2020b+**

No se usan `classdef` folders. Todos los archivos son `.m` individuales.

## Tests

```bash
# Octave
octave --no-gui --eval "run('tests/test_all.m')"

# MATLAB
matlab -batch "run('tests/test_all.m')"
```

## Referencias

- Kogelnik, H., & Li, T. (1966). Laser beams and resonators. *Applied Optics*.
- Siegman, A. E. (1986). *Lasers*. University Science Books.

## Migracion Legacy

- Plan incremental (Strangler): `docs/migration/LEGACY_MIGRATION_PLAN.md`
