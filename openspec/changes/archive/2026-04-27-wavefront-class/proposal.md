# Proposal: Wavefront Class

## Intent

Introduce a `Wavefront` class for analyzing optical wavefronts from beam fields. This enables Zernike decomposition, wavefront error analysis, and Strehl ratio calculation — common tasks in optical system characterization and adaptive optics.

## Scope

### In Scope
- Create `Wavefront` class in `src/visualization/` (not a package yet)
- Store complex field data `[Ny x Nx]` and metadata (wavelength, grid)
- Compute intensity and phase from complex field
- Calculate wavefront gradient (for deformable mirror input)
- Zernike polynomial fitting (standard 36-term set)
- Wavefront error metrics: RMS, PV (peak-to-valley)
- Strehl ratio estimation from wavefront error
- Static visualization methods for wavefront and Zernike coefficients

### Out of Scope
- Adaptive optics closed-loop control
- Temporal wavefront analysis
- Integration with DeformableMirror class (future work)
- Propagation through optical systems

## Capabilities

### New Capabilities
- `wavefront-extraction`: Extract wavefront (phase) from complex field data
- `zernike-fitting`: Fit Zernike polynomials to wavefront error
- `wavefront-metrics`: Compute RMS, PV, Strehl ratio
- `wavefront-visualization`: Plot wavefront maps and Zernike bar charts

## Approach

1. **Input**: Receive complex field `[Ny x Nx]` from existing beam's `opticalField()`
2. **Storage**: Store field, wavelength, grid spacing, reference radius
3. **Analysis**: Phase extraction → unwrapping (if needed) → Zernike fit → metrics
4. **Output**: Metrics struct, Zernike coefficients vector, residual map
5. **Visualization**: 2D phase map, Zernike coefficient bar chart, intensity overlay

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `src/visualization/` | New | New `Wavefront.m` class |
| `docs/ARCHITECTURE.md` | Modified | Add Wavefront section |
| `examples/canonical/` | Modified | Add `ExampleWavefront.m` |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Phase unwrapping edge cases | Medium | Start with analytical beams (known phase) for testing |
| Zernike normalization | Low | Use standard Noll indexing and normalization |
| Grid coordinate handling | Medium | Accept GridUtils input directly |

## Rollback Plan

1. Revert PR, delete `src/visualization/Wavefront.m`
2. Remove `ExampleWavefront.m`
3. Revert ARCHITECTURE.md changes

## Dependencies

- Requires: `GridUtils` for grid metadata
- Works with: any `ParaxialBeam` subclass

## Success Criteria

- [ ] `Wavefront(field, lambda, grid)` constructor works
- [ ] `getPhase()` returns unwrapped phase `[Ny x Nx]`
- [ ] `fitZernike(nTerms)` returns coefficient vector
- [ ] `computeRMS()` and `computePV()` return scalar values
- [ ] `computeStrehl()` returns ratio 0-1
- [ ] `plotWavefront()` displays phase map
- [ ] Test with GaussianBeam field: known analytical solution for verification
- [ ] Tests pass in Octave and MATLAB CI