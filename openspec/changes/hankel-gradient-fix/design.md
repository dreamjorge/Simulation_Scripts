# Design: HankelLaguerre Gradient Explosion Fix

## Technical Approach

**Chosen: Option C — Hybrid approach with polar-coordinate phase gradient for HankelLaguerre beams with l≠0.**

The root cause: HankelLaguerre beams carry phase `exp(i·l·θ)`. In Cartesian coordinates, central-difference evaluation `u(x±δ, y)` lands on different θ values, causing the gradient to sample across the 2π branch cut. This produces explosive, asymmetric slopes at θ=0 and θ=π.

The fix: Compute phase gradient directly in polar coordinates (r, θ), where the vortex structure is natural and the θ-derivative is periodic (no branch cuts). Transform back to Cartesian via:

```
∇φ = ∂φ/∂r · (x̂/r, ŷ/r) + ∂φ/∂θ · (-ŷ/r², x̂/r²)
```

Non-singular beams (Gaussian, Hermite, l=0 HankelLaguerre) continue using the existing Cartesian complex gradient — no performance penalty.

## Architecture Decisions

### Decision: Polar gradient transformation is well-defined at vortex cores

**Choice**: Use `∇φ_cart = ∂φ/∂r · (x/r, y/r) + ∂φ/∂θ · (-y/r², x/r²)` everywhere, including near r=0.

**Alternatives considered**:
- Option B (delta reduction): Impractical — requires knowing phase branch structure a priori to avoid landing on different branches.
- Direct Cartesian fallback near axis: Same branch-cut problem, just at smaller scale.

**Rationale**: The polar-to-Cartesian gradient transformation is mathematically exact. At r→0 with l≠0, the beam amplitude |u|→0 (vortex), so Im{u*·∇u}/|u|² is regularized by the existing ε=10⁻¹² term. The physical ray direction is still well-defined — it's purely azimuthal (tangential to the vortex).

### Decision: Detect vortex singularity via beam type, not field evaluation

**Choice**: `beamHasVortex(beam) = isa(beam,'HankelLaguerre') && beam.l ~= 0`

**Alternatives considered**:
- Monitor field amplitude near evaluation points: Too expensive, triggers falsely for weak sidelobes.
- Check if |u|² < threshold at evaluation points: ε=10⁻¹² regularization already handles this — adding another detection layer adds complexity without benefit.

**Rationale**: The vortex exists only for HankelLaguerre with l≠0. This is a compile-time property we know from the beam's constructor. No runtime detection needed.

### Decision: New method on `RayTracer` for polar-coordinate gradient

**Choice**: Add static method `calculatePhaseGradientPolar(beam, x, y, z)` and call it from `HankelRayTracer.calculateSlopes`.

**Alternatives considered**:
- Modify `calculatePhaseGradientComplex` directly: Would require passing beam type info, violating single-responsibility.
- Create a separate utility class: Overkill — polar/Cartesian gradient is a beam property, not a general utility.

**Rationale**: `HankelRayTracer.calculateSlopes` already dispatches per ht-type. Adding a polar-gradient branch there keeps the vortex handling localized to the Hankel tracer, while `RayTracer.calculateSlopes` remains unchanged for non-Hankel beams.

## Data Flow

```
HankelRayTracer.propagate
  └─> HankelRayTracer.calculateSlopes(beam, x, y, z, ht)
        ├─ For non-vortex beams: RayTracer.calculatePhaseGradientComplex (unchanged)
        └─ For vortex beams (l≠0): RayTracer.calculatePhaseGradientPolar
              ├─ cart2pol(x, y) → (r, θ)
              ├─ Central diff in (r, θ): ∂u/∂r, ∂u/∂θ
              ├─ Phase gradient: Im{u*·∂u/∂r}/|u|², Im{u*·∂u/∂θ}/|u|²
              └─ Transform: sx = (1/k)·[dφ/dr · x/r - dφ/dθ · y/r²]
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/propagation/rays/RayTracer.m` | Modify | Add `calculatePhaseGradientPolar` static method; add `beamHasVortex` static helper |
| `src/propagation/rays/HankelRayTracer.m` | Modify | Dispatch to polar gradient when `beamHasVortex(beam)` is true |
| `tests/edge_cases/test_HankelLaguerre_gradient.m` | Create | Regression test: verify gradient is finite and physical at θ=0, π, and off-axis |

## Algorithm: `calculatePhaseGradientPolar`

```
Input: beam, x, y (matrices), z
Output: sx, sy (slopes dx/dz, dy/dz)

1. Convert to polar: [TH, R] = cart2pol(x, y)

2. Resolve delta_r = max(λ, |x|·10⁻⁴, |y|·10⁻⁴, w₀·10⁻⁴)  [same as resolveDelta]
   delta_theta = delta_r / R  (arc-length equivalent: r·δθ = δr → δθ = δr/r)

3. Evaluate field at 5 points:
   u0   = beam.opticalField(x, y, z)
   u_rp = beam.opticalField(perturbed_in_r_positive_direction, ...)
   u_rm = beam.opticalField(perturbed_in_r_negative_direction, ...)
   u_tp = beam.opticalField(perturbed_in_theta_positive_direction, ...)
   u_tm = beam.opticalField(perturbed_in_theta_negative_direction, ...)

4. Central differences:
   dudr = (u_rp - u_rm) / (2·delta_r)
   dudt = (u_tp - u_tm) / (2·delta_theta)

5. Phase gradients via Im{u*·∇u}/|u|²:
   dphidr = imag(conj(u0)·dudr) / (|u0|² + ε)
   dphidt = imag(conj(u0)·dudt) / (|u0|² + ε)

6. Transform to Cartesian:
   sx = (1/k)·[dphidr·x/R - dphidt·y/R²]
   sy = (1/k)·[dphidr·y/R + dphidt·x/R²]
```

**Note on delta_theta**: When R is small, delta_theta can become large (division by small number). Cap delta_theta at π/4 to avoid sampling more than a quadrant per step.

## Interfaces

### New method: `RayTracer.calculatePhaseGradientPolar`

```matlab
function [sx, sy] = calculatePhaseGradientPolar(beam, x, y, z)
    % POLAR PHASE GRADIENT — For beams with azimuthal phase singularities.
    %
    % For HankelLaguerre beams with l≠0, the phase has exp(i·l·θ) structure.
    % Central difference in Cartesian (x±δ, y±δ) crosses θ boundaries,
    % producing spurious gradients at θ=0 and θ=π.
    %
    % This method computes ∂φ/∂r and ∂φ/∂θ directly in polar coords,
    % where the vortex is natural and θ-derivative is periodic.
    %
    % Gradient transform: ∇φ = ∂φ/∂r · (r̂) + ∂φ/∂θ · (θ̂/r)
    %                   = ∂φ/∂r · (x/r, y/r) + ∂φ/∂θ · (-y/r², x/r²)
```

### New helper: `RayTracer.beamHasVortex`

```matlab
function [hasVortex] = beamHasVortex(beam)
    % Returns true for HankelLaguerre beams with l ≠ 0.
    hasVortex = isa(beam, 'HankelLaguerre') && (beam.l ~= 0);
end
```

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Unit | `calculatePhaseGradientPolar` for l=0 vs l≠0 | Verify finite gradients at θ=0, π vs explosive gradients before fix |
| Unit | `beamHasVortex` | Gaussian→false, l=0→false, l=8→true |
| Unit | Polar→Cartesian transform at r=1e-6 (near-vortex) | Verify no NaN/Inf |
| Integration | `HankelRayTracer.propagate` with l=8, R_obs=60µm | Verify sy_rk4 ≈ 0 at θ=0 (not 33,300 rad/m) |
| Regression | `HankelRayTracer` with Hermite, l=0 Laguerre | Verify existing correct behavior unchanged |

### Critical Test Cases

1. **HankelLaguerre l=8, p=0, R=60µm, θ=0**: sy_rk4 must be finite and physically plausible (|sy| < 1000 rad/m, not 33,300)
2. **Same beam, θ=π/4**: sy_rk4 ≈ 0 (already works, must not regress)
3. **Same beam, θ=π**: Same behavior as θ=0
4. **GaussianBeam at same geometry**: Unchanged correct behavior
5. **HankelLaguerre l=0**: Same as Gaussian (no vortex) — verify no regression

## Open Questions

- [ ] **Polar delta cap**: Is π/4 (90°) cap on delta_theta sufficient? For very high l (l≥20), even a quarter-rotation may cross branch cuts. Should cap be π/(2·l) to ensure we never sample more than 1/l of a full vortex rotation?
- [ ] **Performance**: The polar method requires 5 `opticalField` calls vs 5 for Cartesian. Is there a way to reuse field evaluations across r and θ perturbations? (Low priority — 5 calls is already minimal)
- [ ] **Hermite beams**: `HankelHermite` does not have a vortex (Cartesian coordinate system). Verify `beamHasVortex` returns false for `HankelHermite` — it should, since `isa(beam,'HankelLaguerre')` will be false.

## Risk Assessment

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Polar gradient still explosive for very high l | Low | Cap delta_theta; add test for l=16, l=32 |
| Division by R² in θ-term causes Inf at r→0 | Medium | ε-regularization in denominator; test at r=1e-9 |
| Breaking existing HankelHermite behavior | Low | `beamHasVortex` returns false for Hermite; existing Cartesian path used |
