# Tasks: raytracing-numerical-improvement

## Phase 1: Core Math — RayTracer Gradient Improvement

- [x] 1.1 Add `calculatePhaseGradientComplex(beam, x, y, z)` private method to `src/propagation/rays/RayTracer.m` — implement `∇φ = Im(u̅∇u) / (|u|² + ε)` with central difference, ε=1e-12
- [x] 1.2 Add `resolveDelta(x, y, w0, lambda)` private static method to `RayTracer.m` — return `max(lambda, abs(x)*1e-4, abs(y)*1e-4, w0*1e-4)`
- [x] 1.3 Replace `calculateSlopes()` body to call `calculatePhaseGradientComplex` instead of forward-diff+unwrap path
- [ ] 1.4 Verify existing `test_RayTracing.m` slopes-at-center test passes with new method

## Phase 2: State Consistency — RayBundle

- [ ] 2.1 Convert `r` and `theta` from stored properties to Dependent in `src/propagation/rays/RayBundle.m`
- [ ] 2.2 Add `get.r()` and `get.theta()` methods computing from `x(:,:,end)` and `y(:,:,end)`
- [ ] 2.3 Remove `r` and `theta` assignment from constructor and `addStep()`
- [ ] 2.4 Verify `tests/modern/test_RayTracing.m` bundle initialization tests still pass

## Phase 3: Hankel Axis-Crossing Fix

- [ ] 3.1 Apply same complex gradient to `HankelRayTracer.calculateSlopes()` — replace forward-diff+unwrap
- [ ] 3.2 Replace axis-crossing determinant `(x0.*y1 - x1.*y0) < 0` with geometric minimum-distance logic in `HankelRayTracer.propagate()`
- [ ] 3.3 Add test in `tests/modern/test_HankelRayTracing.m` verifying flip only on real axis crossing, not on orientation change

## Phase 4: Physical Accuracy Tests

- [ ] 4.1 Add test in `tests/modern/test_RayTracing.m` comparing numerical gradient vs analytical Gaussian gradient `dφ/dx = k·x/R(z)` with tolerance 1e-6
- [ ] 4.2 Add test verifying Euler vs RK4 error decreases as dz refines (convergence order check)
- [ ] 4.3 Add test verifying radial symmetry: `bundle.r` matches `sqrt(x²+y²)` after propagation
- [ ] 4.4 Update `tests/edge_cases/test_RayTracing_extreme.m` — remove or rewrite tests that set `bundle.sx/sy` expecting them to affect propagation (they don't)

## Phase 5: Integration Verification

- [ ] 5.1 Run full `test_RayTracing.m` suite — all tests pass
- [ ] 5.2 Run full `test_HankelRayTracing.m` suite — all tests pass
- [ ] 5.3 Run full `test_RayTracing_extreme.m` suite — no regressions
- [ ] 5.4 Verify gradient accuracy improvement: numerical gradient should deviate < 1e-4 rad/m from analytical for Gaussian at z=0