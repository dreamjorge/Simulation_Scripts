# Simulation_Scripts Test Suite

This directory contains tests for the ParaxialBeams simulation library.

## Running Tests

```bash
cd <repo-root>
octave --no-gui --eval "run('tests/test_all.m')"
```

With full MATLAB installed:
```bash
matlab -batch "run('tests/test_all.m')"
```

Note: MATLAB Runtime alone is not enough for the `matlab -batch` command above.

## Critical Coverage Gates

- `z = 0` no produce NaN/Inf invalidos en tests individuales de beams (`test_GaussianBeam.m`, `test_HermiteBeam.m`, `test_LaguerreBeam.m`, `test_ElegantHermiteBeam.m`, `test_ElegantLaguerreBeam.m`)
- modo cero Hermite/Laguerre mantiene equivalencia razonable con Gaussian en `test_HermiteBeam.m` y `test_LaguerreBeam.m`
- ray tracing cilindrico sigue estable en `test_RayTracing.m`
- `tests/test_all.m` corre via `portable_runner()`
- compatibilidad MATLAB requiere ejecutar la misma suite en un entorno con MATLAB completo

## Current Runner Scope

- `portable_runner()` ejecuta hoy `test_PhysicalConstants.m` y `test_RayTracing.m`
- los tests de beams existen y cubren checks criticos, pero todavia no forman parte del runner portable principal
- `runAllTests.m` y otros helpers quedan como scripts auxiliares, no como entrypoint auditado principal

## Test Coverage

### PhysicalConstants
- Wave number calculation
- Rayleigh distance
- Waist at z position
- Radius of curvature
- Gouy phase

### GridUtils
- 2D grid creation
- Frequency grid generation
- Static meshgrid methods

### FFTUtils
- FFT roundtrip
- Transfer function
- Propagation

### GaussianParameters
- Rayleigh distance, waist, wave number
- Vector z-coordinate support
- Gouy phase
- Divergence angle

### Hermite and Laguerre Models
- HermiteParameters, LaguerreParameters
- GaussianBeam, HermiteBeam, LaguerreBeam
- ElegantHermiteBeam, ElegantLaguerreBeam

### AnalysisUtils
- **gradientRZ**: Calculate ray slope in r-z plane
  - Interior point calculation with 10x10 field slices
  - Handles zero gradient (uniform fields) gracefully
  
- **gradientXYZ**: Calculate local gradients in 3D
  - Returns three gradient values (mzx, mzy, mxy)
  - Bounds checking with 64x64 field slices

### ElegantHermiteParameters
- Stores n and m mode indices
- Alpha complex parameter at z>0
- Alpha at waist matches k/(2*zr)
- Default indices (n=0, m=0)

### ElegantLaguerreParameters
- Stores l (topological charge) and p (radial index)
- Alpha matches ElegantHermiteParameters formula
- Default indices (l=0, p=0)

## Test Patterns

Tests follow consistent patterns:
- Use `fprintf` for PASS/FAIL output
- Increment `passed` and `failed` counters
- Tolerance: 1e-5 for physical values, 1e-10 for phase values
- Octave/MATLAB compatible classdef syntax
