# Delta for ray-slope-calculation

## MODIFIED Requirements

### Requirement: Numerical Phase Gradient Computation

The system SHALL compute spatial phase gradients from the optical field for ray slope determination. The gradient MUST be calculated using a method that is robust to phase singularities and low-amplitude regions.

(Previously: forward difference with fixed delta=1e-7 and unwrap(angle(field)))

#### Scenario: Gradient calculation at beam center

- GIVEN a GaussianBeam at waist (z=0) with waist `w0` and wavelength `λ`
- WHEN `RayTracer.calculateSlopes(beam, x, y, z)` is called at a point where `|field| > ε`
- THEN the returned slopes `(sx, sy)` MUST be numerically stable and finite
- AND the gradient magnitude MUST be consistent with paraxial propagation from phase curvature

#### Scenario: Gradient calculation near phase singularities

- GIVEN a beam with phase discontinuities or vortices (phase undefined at origin)
- WHEN `calculateSlopes` is called at a point near the singularity where `|field| > ε`
- THEN the method MUST NOT rely on `unwrap(angle(field))` for gradient computation
- AND the result MUST be consistent with the physical field structure

#### Scenario: Gradient calculation in low-amplitude regions

- GIVEN a point where `|field| < ε` (far from beam center or in amplitude null)
- WHEN `calculateSlopes` is called
- THEN the system MUST apply regularization to avoid division by zero
- AND return a fallback gradient based on local field behavior

#### Scenario: Delta resolution adapts to spatial scale

- GIVEN a beam with characteristic scale `w0` and wavelength `λ`
- WHEN `calculateSlopes` is called at position `(x, y)`
- THEN the perturbation delta MUST be computed as `resolveDelta(x, y, w0, λ)` returning `max(λ, |x|*1e-4, |y|*1e-4)`
- AND NOT be hardcoded to a fixed value

#### Scenario: Slope calculation uses complex field directly

- GIVEN an optical field `u(x,y,z)` with amplitude and phase
- WHEN the system computes phase gradients
- THEN it SHOULD compute `∇φ = Im(u̅ ∇u) / |u|²` where `u̅` is the complex conjugate of `u`
- OR fall back to central difference if the complex method fails

#### Scenario: Propagation integrates slopes correctly

- GIVEN a RayBundle with initial coordinates `(x₀, y₀, z₀)` and computed slopes `(sx, sy)`
- WHEN `RayTracer.propagate` is called with integration method Euler or RK4
- THEN the ray positions at subsequent z steps MUST follow the integrated trajectories
- AND the slope computation at each step MUST use the local field at the current ray position