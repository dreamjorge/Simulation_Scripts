# Delta for complex-phase-gradient

## ADDED Requirements

### Requirement: Complex Field Phase Gradient

The system SHALL compute phase gradients directly from the complex optical field without phase unwrapping, using the identity `∇φ = Im(u̅ ∇u) / |u|²`.

The complex field method MUST:
- Evaluate the field and its spatial derivatives via central difference
- Use `u̅ * ∂u/∂x` and `u̅ * ∂u/∂y` (conjugate times gradient)
- Normalize by `|u|² + ε` where `ε` is a small regularization constant
- Fall back to central difference phase gradient if direct method fails

#### Scenario: Complex gradient at arbitrary position

- GIVEN a beam implementing `opticalField(x, y, z)` and any beam parameters
- WHEN `calculatePhaseGradientComplex(beam, x, y, z)` is invoked
- THEN the method MUST return `(sx, sy)` as the normalized imaginary part of conjugate times field gradient
- AND the result MUST have units of radians per meter, matching `dx/dz` for paraxial rays

#### Scenario: Regularization prevents division by zero

- GIVEN a field amplitude `|field| < ε` where `ε = 1e-12`
- WHEN the complex gradient is computed
- THEN the denominator MUST use `|field|² + ε` instead of pure `|field|²`
- AND return a numerically stable gradient even in amplitude nulls

#### Scenario: Central difference fallback for spatial derivatives

- GIVEN the need to compute `∂u/∂x` and `∂u/∂y` for the complex gradient
- WHEN `resolveDelta(x, y, w0, λ)` returns a local delta
- THEN the derivative MUST be computed via central difference: `∂u/∂x ≈ (u(x+δ) - u(x-δ)) / (2δ)`
- AND the same for y direction with symmetric delta

#### Scenario: Resolution adapts to geometry

- GIVEN a ray at position `(x, y)` with local spatial scale determined by beam waist `w0`
- WHEN `resolveDelta` is called
- THEN the returned delta MUST be `max(λ, |x|*1e-4, |y|*1e-4, w0*1e-4)`
- AND the delta MUST be large enough to exceed numerical noise but small enough to capture field curvature

#### Scenario: Integration method determines slope aggregation

- GIVEN the RK4 integration method requires slope evaluations at intermediate points
- WHEN `calculateSlopes` is called multiple times with offset positions during RK4 stages
- THEN each evaluation MUST use the same delta resolution logic for consistency
- AND the final slope for the step MUST be the proper RK4 weighted average of stage slopes

### Requirement: RayBundle Polar Coordinates as Dependent

The system SHALL compute polar coordinates `(r, θ)` from Cartesian `(x, y)` on every access, ensuring consistency at all times.

#### Scenario: r and theta update after addStep

- GIVEN a RayBundle that has completed several `addStep()` calls
- WHEN `bundle.r` or `bundle.theta` is accessed
- THEN the values MUST be computed from the current `bundle.x` and `bundle.y` matrices
- AND MUST reflect the latest ray positions, not stale values from initialization

#### Scenario: addStep does not store r and theta

- GIVEN a RayBundle with `r` and `theta` declared as Dependent
- WHEN `addStep(x, y, z, sx, sy)` is called
- THEN the implementation MUST NOT store `r` or `theta` in the step
- AND subsequent access MUST recompute from x, y at the latest z-slice

### Requirement: Axis Crossing Detection via Minimum Distance

The system SHALL determine optical axis crossings using geometric distance from ray segment to origin, not orientation sign changes.

#### Scenario: Real axis crossing flips Hankel type

- GIVEN a HankelLaguerre ray at position `(x₀, y₀)` approaching the origin
- WHEN the ray segment `(x₀,y₀) → (x₁,y₁)` passes within minimum distance `d < threshold` of origin
- AND the orthogonal projection of origin onto the segment falls within `[0,1]` of the segment parameter
- THEN the Hankel type MUST flip from `H^(2)` to `H^(1)` when crossing occurs

#### Scenario: Near-miss does not flip type

- GIVEN a ray segment that passes near but not through the origin
- WHEN the minimum distance from segment to origin exceeds the threshold
- THEN the Hankel type MUST NOT flip regardless of orientation sign of `(x₀*y₁ - x₁*y₀)`

#### Scenario: Grazing angle crossing detection

- GIVEN a ray with small angle relative to optical axis that closely approaches origin
- WHEN the segment passes within `max(|x₀|, |y₀|) * 1e-3` of origin
- THEN crossing MUST be detected even if neither endpoint is at origin
- AND the flip MUST occur at the point of closest approach

#### Scenario: No spurious flip for rays not aimed at origin

- GIVEN a ray trajectory that curves away from origin
- WHEN the segment endpoint is farther from origin than the start point
- THEN even if the orientation sign changes, the Hankel type MUST NOT flip