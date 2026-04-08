# Simulation_Scripts Test Suite

This directory contains tests for the ParaxialBeams simulation library.

## Running Tests

```bash
cd /root/Simulation_Scripts
octave --no-gui --eval "run('tests/test_all.m')"
```

Or in MATLAB:
```matlab
run('tests/test_all.m')
```

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
