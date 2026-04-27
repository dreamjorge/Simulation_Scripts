# Tasks: Wavefront Class

## Phase 1: Zernike Utilities

- [x] 1.1 Create `src/visualization/ZernikeUtils.m` with static `zernike(n, rho, theta)` method
- [x] 1.2 Implement Noll Zernike formulas for n=1..36
- [x] 1.3 Add `zernikeName(n)` static method returning names: 'Piston', 'Tilt X', 'Defocus', etc.
- [x] 1.4 Add `zernikeMatrix(rho, theta, nTerms)` static method building design matrix

## Phase 2: Wavefront Core Class

- [x] 2.1 Create `src/visualization/Wavefront.m` classdef with properties: field, wavelength, dx, dy, Dx, Dy
- [x] 2.2 Implement constructor: `Wavefront(E, lambda)` and `Wavefront(E, lambda, grid)`
- [x] 2.3 Implement `getField()` returning complex field
- [x] 2.4 Implement `getIntensity()` returning `abs(E).^2`
- [x] 2.5 Implement `getPhase()` returning `angle(E)` (wrapped)

## Phase 3: Grid Helpers

- [x] 3.1 Add `gridPolar()` method returning rho, theta from stored grid info
- [x] 3.2 Add `gridCartesian()` method returning X, Y grids

## Phase 4: Zernike Fitting

- [x] 4.1 Implement `fitZernike(nTerms)` using `pinv(ZernikeMatrix)` least-squares
- [x] 4.2 Implement `reconstructZernike(coeffs)` to rebuild phase from coefficients
- [x] 4.3 Implement `zernikeResidual(nTerms)` returning residual error

## Phase 5: Metrics

- [x] 5.1 Implement `computeRMS()` returning scalar RMS of phase
- [x] 5.2 Implement `computePV()` returning `max(phi) - min(phi)`
- [x] 5.3 Implement `computeStrehl()` using Maréchal approximation: `exp(-sigma^2)`
- [x] 5.4 Implement `getMetrics()` returning struct with rms, pv, strehl

## Phase 6: Visualization

- [x] 6.1 Implement `plotWavefront()` showing 2D phase map with colorbar
- [x] 6.2 Implement `plotIntensity()` showing 2D intensity map
- [x] 6.3 Implement `plotZernikeCoeffs(coeffs)` bar chart with Zernike names
- [x] 6.4 Implement `plotPhaseSlice(plane, idx)` 1D cross-section

## Phase 7: Tests

- [x] 7.1 Create `tests/modern/test_Wavefront.m` with basic constructor tests
- [x] 7.2 Test Zernike Z1=1 (piston), Z2∝cos (tilt X), Z4=sqrt(3)*(2ρ²-1) (defocus)
- [x] 7.3 Test: Gaussian beam at waist → only piston term ≠ 0 (all others ~0)
- [x] 7.4 Test: fitZernike → reconstructZernike round-trip residual < 1e-10
- [x] 7.5 Test: RMS/PV calculations against known values
- [x] 7.6 Test: Strehl calculation for sigma = 0.1 rad → ~0.990

## Phase 8: Documentation & Examples

- [x] 8.1 Create `examples/canonical/ExampleWavefront.m` demonstrating all features
- [x] 8.2 Update `docs/ARCHITECTURE.md` with Wavefront section
- [x] 8.3 Register `tests/modern/test_Wavefront.m` in `tests/portable_runner.m`