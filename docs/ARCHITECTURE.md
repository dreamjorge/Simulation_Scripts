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

## Phase Convention

- The beam implementations use the phasor convention `exp(-i*k*z)`.
- Gaussian carrier phase terms are `exp(-i*psi(z))` for Gouy and `exp(+i*k*r^2/(2R(z)))` for curvature.
- Higher-order modal terms use the same sign convention as the carrier Gouy term (`exp(-i*phi_mode)`), where:
  - Hermite: `phi_mode = (n+m)*psi(z)`
  - Laguerre/Hankel: `phi_mode = (|l|+2p)*psi(z)`

## Data Flow

```
┌──────────────────────────────────────────────────────────────────────┐
│                         User Script                                   │
│   examples/canonical/MainGauss_refactored.m, etc.                    │
└──────────────────────────────┬───────────────────────────────────────┘
                               │ setpaths / addpath src + ParaxialBeams
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
├── src/
│   ├── beams/
│   │   ├── ParaxialBeam.m        ← Abstract base (contract)
│   │   ├── GaussianBeam.m
│   │   ├── HermiteBeam.m
│   │   ├── LaguerreBeam.m
│   │   ├── ElegantHermiteBeam.m
│   │   ├── ElegantLaguerreBeam.m
│   │   └── HankelLaguerre.m
│   ├── parameters/
│   │   ├── GaussianParameters.m
│   │   ├── HermiteParameters.m
│   │   ├── LaguerreParameters.m
│   │   ├── ElegantHermiteParameters.m
│   │   └── ElegantLaguerreParameters.m
│   ├── propagation/
│   │   ├── field/
│   │   │   ├── IPropagator.m         ← Strategy interface
│   │   │   ├── FFTPropagator.m
│   │   │   └── AnalyticPropagator.m
│   │   └── rays/
│   │       ├── RayTracePropagator.m
│   │       ├── OpticalRay.m
│   │       ├── CylindricalRay.m
│   │       ├── RayBundle.m
│   │       └── RayTracer.m
│   └── visualization/
│       └── VisualizationUtils.m
├── ParaxialBeams/
│   ├── BeamFactory.m         ← Factory
│   ├── PhysicalConstants.m
│   ├── GridUtils.m
│   ├── FFTUtils.m
│   ├── AnalysisUtils.m
│   ├── PolynomialUtils.m
│   └── Addons/               ← Helper functions
├── examples/                  ← Usage examples
│   ├── canonical/
│   │   ├── MainGauss_refactored.m
│   │   ├── MainMultiMode.m
│   │   └── ExampleRayTracing.m
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

## Referencias Bibliograficas por Implementacion

### GaussianBeam (`src/beams/GaussianBeam.m`)

- Kogelnik, H., & Li, T. (1966). *Laser Beams and Resonators*. Applied Optics, 5(10), 1550-1567.
  - Base para el formalismo paraxial: waist `w(z)`, radio de curvatura `R(z)`, fase de Gouy `psi(z)`.
- Siegman, A. E. (1986). *Lasers*. University Science Books.
  - Referencia clasica para notacion y convenciones de fase de haces gaussianos.

### HermiteBeam / LaguerreBeam (`src/beams/HermiteBeam.m`, `src/beams/LaguerreBeam.m`)

- Saleh, B. E. A., & Teich, M. C. (2007). *Fundamentals of Photonics* (2nd ed.). Wiley.
  - Derivacion de modos transversales HG/LG y estructura modal en coordenadas cartesianas/polares.
- Allen, L., Beijersbergen, M. W., Spreeuw, R. J. C., & Woerdman, J. P. (1992).
  *Orbital angular momentum of light and the transformation of Laguerre-Gaussian laser modes*.
  Physical Review A, 45(11), 8185-8189.
  - Marco para indice azimutal `l`, orden radial `p` y fase azimutal `exp(i*l*theta)`.

### ElegantHermiteBeam / ElegantLaguerreBeam (`src/beams/ElegantHermiteBeam.m`, `src/beams/ElegantLaguerreBeam.m`)

- Siegman, A. E. (1996). *Defining and measuring laser beam quality*. JOSA A, 13(5), 952-964.
  - Uso del parametro complejo `q(z)` y formulacion "elegant" con argumentos complejos.
- Siegman, A. E. (1986). *Lasers*. University Science Books.
  - Contexto teorico para el acople amplitud-fase en familias elegant.

### HankelLaguerre (`src/beams/HankelLaguerre.m`)

- Kotlyar, V. V., Kovalev, A. A., & Porfirev, A. P. (2012). *Hankel beams and their properties*.
  Optics Letters / JOSA A literature on Hankel-type structured beams.
  - Base conceptual para la combinacion `H^(1,2) = LG +/- i*XLG` usada en la implementacion.

### Soporte MATLAB/Octave (semantica numerica)

- MathWorks. `cart2pol` documentation.
  - Verifica orden de salida `[theta, rho]`, usado al mapear `(X,Y)->(theta,r)`.
- MathWorks. Symbolic Math Toolbox `laguerreL` documentation.
  - Convencion de polinomios de Laguerre asociados consistente con `L_p^{|l|}`.

### Nota de Convenciones

- Distintas fuentes usan convenciones de fasor distintas (`exp(+i*k*z)` vs `exp(-i*k*z)`).
- Este repositorio fija `exp(-i*k*z)` y mantiene coherencia interna en:
  - carrier gaussiano,
  - fase de curvatura,
  - fase de Gouy modal (`exp(-i*phi_mode)`).
