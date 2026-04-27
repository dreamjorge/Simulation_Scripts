# Design: Wavefront Class

## Technical Approach

Create `Wavefront` as a standalone class in `src/visualization/` that wraps complex field data and provides wavefront analysis methods. It follows the existing class patterns (static helpers, static computation delegation).

## Architecture Decisions

### Decision: Class Location

**Choice**: Create in `src/visualization/Wavefront.m`

**Alternatives**: New `src/wavefront/` folder
**Rationale**: `Wavefront` is primarily a visualization/analysis tool; fits with `VisualizationUtils`.

### Decision: Zernike Implementation

**Choice**: Use Noll indexing with standard definitions; compute via explicit formulas

```matlab
methods (Static)
    function Z = zernike(n, rho, theta)
        % Noll standard Zernike polynomials
        switch n
            case 1, Z = ones(size(rho));                    % Piston
            case 2, Z = 2*rho.*cos(theta);                  % Tilt X
            case 3, Z = 2*rho.*sin(theta);                  % Tilt Y
            case 4, Z = sqrt(3)*(2*rho.^2 - 1);            % Defocus
            % ... up to 36
        end
    end
end
```

**Alternatives**: Use symbolic toolbox, external library
**Rationale**: Explicit formulas are self-contained, no dependencies, fast.

### Decision: Least-Squares Fitting

**Choice**: Use pseudo-inverse for Zernike coefficient fitting

```matlab
function coeffs = fitZernike(obj, nTerms)
    [rho, theta] = obj.gridPolar();
    Z = obj.buildZernikeMatrix(rho, theta, nTerms);
    coeffs = pinv(Z) * obj.phase(:);  % Least-squares
end
```

**Rationale**: `pinv` handles ill-conditioned matrices gracefully.

### Decision: Phase Storage

**Choice**: Store wrapped phase (`angle(E)`); unwrap on demand

**Rationale**: Preserves original data; unwrapping is lossy.

## Data Flow

```
ParaxialBeam.opticalField(X, Y, z)
        │
        ▼ (complex field)
   Wavefront(field, lambda, grid)
        │
        ├── getPhase() ────────────────────► wrap/unwrap
        ├── fitZernike(n) ─────────────────► Zernike coefficients
        ├── computeRMS() ──────────────────► scalar metrics
        ├── computeStrehl() ───────────────► Strehl ratio
        └── plotWavefront() ───────────────► visualization
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/visualization/Wavefront.m` | Create | Main Wavefront class |
| `src/visualization/ZernikeUtils.m` | Create | Static Zernike computation helpers |
| `docs/ARCHITECTURE.md` | Modify | Add Wavefront section |
| `examples/canonical/ExampleWavefront.m` | Create | New canonical example |

## Interfaces / Contracts

### Wavefront Constructor

```matlab
wf = Wavefront(E, lambda)           % Minimal
wf = Wavefront(E, lambda, grid)     % With GridUtils
wf = Wavefront(E, lambda, dx, dy)   % With explicit spacing
```

### Key Methods

```matlab
phi = wf.getPhase()                    % Wrapped phase
I = wf.getIntensity()                  % |E|^2
coeffs = wf.fitZernike(36)            % Fit 36 terms
rms = wf.computeRMS()                  % RMS wavefront error
pv = wf.computePV()                    % Peak-to-valley
strehl = wf.computeStrehl()            % Strehl ratio
wf.plotWavefront()                     % Phase map visualization
```

### ZernikeUtils Static Helpers

```matlab
ZernikeUtils.zernike(n, rho, theta)      % Compute Z_n
ZernikeUtils.zernikeName(n)               % 'Defocus', etc.
ZernikeUtils.zernikeMatrix(rho, theta, N) % Build design matrix
```

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Unit | Zernike formulas vs known values | Test Z1=1, Z2∝cos, etc. |
| Unit | RMS/PV calculations | GaussianBeam has analytic solution |
| Unit | Strehl formula | Verify exp(-(2πσ/λ)²) |
| Integration | Round-trip: fit → reconstruct | RMS residual < 1e-10 |
| Integration | GaussianBeam → Wavefront | Known analytical wavefront |

### Verification Test: GaussianBeam

Gaussian beam at waist (z=0) has planar wavefront → only piston term ≠ 0.

```matlab
beam = BeamFactory.create('gaussian', 100e-6, 632.8e-9);
E = beam.opticalField(X, Y, 0);
wf = Wavefront(E, lambda, grid);
coeffs = wf.fitZernike(36);
% coeffs(1) ≈ mean phase, all others ≈ 0
```

## Open Questions

- [ ] Unwrapping algorithm: use IED (Itoh) or Goldstein? Recommendation: **Goldstein** for simplicity
- [ ] Should Wavefront support multiple wavelengths per instance? Recommendation: **No** — one wavelength per instance
- [ ] Zernike normalization: Noll standard vs. Fringe? Recommendation: **Noll** (standard in optics)