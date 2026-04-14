# Simulation_Scripts Architecture

## Overview

MATLAB/Octave library for paraxial optical beam propagation simulation. Supports Gaussian, Hermite-Gaussian, Laguerre-Gaussian, and Elegant beam modes with multiple propagation methods.

## Class Hierarchy

```
ParaxialBeam (abstract base)
├── GaussianBeam
├── HermiteBeam
├── LaguerreBeam
├── ElegantHermiteBeam
├── ElegantLaguerreBeam
└── HankelLaguerre
```

## Design Patterns

### 1. Strategy Pattern — Propagators

All propagation algorithms implement the `IPropagator` interface, enabling runtime selection of propagation method.

```
IPropagator (interface)
├── FFTPropagator       — Angular spectrum via FFT
├── AnalyticPropagator  — Direct beam.opticalField evaluation
└── RayTracePropagator  — Phase-gradient ray tracing
```

**Why**: Users can swap propagation methods without changing beam code.

**Usage**:
```matlab
beam = GaussianBeam(100e-6, 632.8e-9);
grid = GridUtils(256, 256, 1e-3, 1e-3);

prop = FFTPropagator(grid, 632.8e-9);
field = prop.propagate(beam, 0.1);
```

### 2. Factory Pattern — BeamFactory

`BeamFactory.create()` instantiates any beam type by name string.

```matlab
% Create by type string
beam = BeamFactory.create('gaussian', 100e-6, 632.8e-9);
beam = BeamFactory.create('hermite', 100e-6, 632.8e-9, 'n', 2, 'm', 1);
beam = BeamFactory.create('laguerre', 100e-6, 632.8e-9, 'l', 1, 'p', 0);
```

**Why**: Decouples calling code from concrete beam classes, simplifies API.

### 3. Parameters Pattern — BeamParameters

Each beam type has a corresponding Parameters class:

```
GaussianParameters
HermiteParameters
LaguerreParameters
ElegantHermiteParameters
ElegantLaguerreParameters
```

Parameters classes provide:
- Waist at given z
- Rayleigh distance
- Gouy phase
- Radius of curvature

**Why**: Separates beam geometry (parameters) from field computation.

## Beam API Contract

Every `ParaxialBeam` subclass MUST implement:

| Method | Returns | Description |
|--------|---------|-------------|
| `opticalField(X, Y, z)` | `[Ny x Nx] complex` | Complex field on Cartesian grid |
| `getParameters(z)` | `GaussianParameters` | Beam params evaluated at z |
| `beamName()` | `char` | String identifier like 'hermite_3_2' |

## Data Flow

```
┌──────────────────────────────────────────────────────────────────────┐
│                         User Script                                   │
│   examples/MainGauss_refactored.m, MainMultiMode.m, etc.             │
└──────────────────────────────┬───────────────────────────────────────┘
                               │ addpath ParaxialBeams
                               ▼
┌──────────────────────────────────────────────────────────────────────┐
│                         BeamFactory                                   │
│        beam = BeamFactory.create(type, w0, lambda, ...)              │
└──────────────────────────────┬───────────────────────────────────────┘
                               │ creates concrete beam
                               ▼
┌──────────────────────────────────────────────────────────────────────┐
│                    ParaxialBeam (abstract)                            │
│   opticalField(X,Y,z) │ getParameters(z) │ beamName()                │
└───────┬───────┬───────┬───────┬───────┬───────┬───────────────────────┘
        │       │       │       │       │       │
   ┌────┴┐ ┌───┴───┐ ┌─┴──┐ ┌──┴──┐ ┌──┴──┐ ┌────┴────┐
   │Gauss│ │Hermit│ │Lag │ │ElegH│ │ElegL│ │HankelL │
   └─────┘ └──────┘ └────┘ └─────┘ └─────┘ └────────┘
        │       │       │       │       │       │
        └───────┴───────┴───────┴───────┴───────┘
                          │ propagate(beam, z)
                          ▼
┌──────────────────────────────────────────────────────────────────────┐
│                      IPropagator (Strategy)                          │
│                                                                       │
│   ┌────────────────┐ ┌─────────────────┐ ┌─────────────────────┐   │
│   │ FFTPropagator  │ │AnalyticPropagat│ │RayTracePropagator   │   │
│   │ angular spec.  │ │ direct formula │ │ phase-gradient RK4 │   │
│   └────────────────┘ └─────────────────┘ └─────────────────────┘   │
└──────────────────────────────────────────────────────────────────────┘
                          │
                          ▼
┌──────────────────────────────────────────────────────────────────────┐
│                        Output                                         │
│   complex field [Ny x Nx]  or  RayBundle (for ray tracing)            │
└──────────────────────────────────────────────────────────────────────┘
```

## Utility Classes

| Class | Purpose |
|-------|---------|
| `PhysicalConstants` | k, zr, waist, R, Gouy phase calculations |
| `GridUtils` | 2D grid, frequency grid, polar grid creation |
| `FFTUtils` | Normalized FFT operations |
| `AnalysisUtils` | Phase gradient to ray slope conversion |
| `PolynomialUtils` | Hermite and Laguerre polynomial evaluation |
| `VisualizationUtils` | Plot helpers |

## Directory Structure

```
Simulation_Scripts/
├── ParaxialBeams/
│   ├── ParaxialBeam.m        ← Abstract base (contract)
│   ├── GaussianBeam.m
│   ├── HermiteBeam.m
│   ├── LaguerreBeam.m
│   ├── ElegantHermiteBeam.m
│   ├── ElegantLaguerreBeam.m
│   ├── HankelLaguerre.m
│   ├── GaussianParameters.m
│   ├── HermiteParameters.m
│   ├── LaguerreParameters.m
│   ├── ElegantHermiteParameters.m
│   ├── ElegantLaguerreParameters.m
│   ├── IPropagator.m         ← Strategy interface
│   ├── FFTPropagator.m
│   ├── AnalyticPropagator.m
│   ├── RayTracePropagator.m
│   ├── BeamFactory.m         ← Factory
│   ├── PhysicalConstants.m
│   ├── GridUtils.m
│   ├── FFTUtils.m
│   ├── AnalysisUtils.m
│   ├── PolynomialUtils.m
│   ├── VisualizationUtils.m
│   ├── OpticalRay.m
│   ├── CylindricalRay.m
│   ├── RayBundle.m
│   ├── RayTracer.m
│   └── Addons/               ← Helper functions
├── examples/                  ← Usage examples
│   ├── MainGauss_refactored.m  %% canonical
│   ├── MainMultiMode.m         %% canonical
│   ├── ExampleRayTracing.m      %% canonical
│   └── ... (many legacy scripts)
├── tests/                     ← Test suite (~380 tests)
├── docs/
│   └── ARCHITECTURE.md
├── README.md
└── plan.md
```

## Dependencies

```
PhysicalConstants (singleton utility)
     │
     └── ParaxialBeam.k, GridUtils.createFreqGrid()

GridUtils ──────────────────→ FFTUtils (frequency grid)

ParaxialBeam ────────────────→ IPropagator (Strategy pattern)
     │                               │
     │                    ┌──────────┼──────────┐
     │                    ▼          ▼          ▼
     │              FFTProp   Analytic   RayTrace
     │
     └── getParameters ──→ GaussianParameters/HermiteParameters/...
```

## Octave/MATLAB Compatibility

All code is written for compatibility with:
- **GNU Octave 11.1.0+**
- **MATLAB R2020b+**

No `classdef` folders (@Folder). All classes in single `.m` files.

## References

- Kogelnik, H., & Li, T. (1966). Laser beams and resonators. *Applied Optics*.
- Siegman, A. E. (1986). *Lasers*. University Science Books.
