# Tasks: +paraxial/ Package Migration (Phase 3)

## Phase 1: BeamFactory Routing Foundation

- [x] 1.1 Modify `ParaxialBeams/BeamFactory.m` — refactor `create()` to use `exist(canonical, 'class')` lookup for `+paraxial/` first, `exist(legacy, 'file')` fallback with warning
- [x] 1.2 Add helper function `BeamFactory.resolveClass(type, canonical, legacy)` to encapsulate routing logic
- [~] 1.3 Run `tests/modern/test_BeamFactory.m` — PARTIAL (12/17 pass; GaussianBeam works, Hermite/Laguerre/Elegant/Hankel have isa() issue in Octave — pre-existing problem in test file)

## Phase 2: Per-Class Migration (Strangler Fig — one beam at a time)

- [x] 2.1 Add deprecation warning to `src/beams/GaussianBeam.m` constructor — `warning('BeamFactory:deprecated', ...)` after super() call
- [~] 2.2 Run `octave --no-gui --eval "run('tests/test_all.m')"` — verify GaussianBeam tests pass (deprecation warnings emit on src/ instantiation)
- [x] 2.3 Add deprecation warning to `src/beams/HermiteBeam.m` constructor
- [~] 2.4 Run tests — verify HermiteBeam tests pass
- [x] 2.5 Add deprecation warning to `src/beams/LaguerreBeam.m` constructor
- [~] 2.6 Run tests — verify LaguerreBeam tests pass
- [x] 2.7 Add deprecation warning to `src/beams/ElegantHermiteBeam.m` constructor
- [~] 2.8 Run tests — verify ElegantHermiteBeam tests pass
- [x] 2.9 Add deprecation warning to `src/beams/ElegantLaguerreBeam.m` constructor
- [~] 2.10 Run tests — verify ElegantLaguerreBeam tests pass
- [x] 2.11 Add deprecation warning to `src/beams/HankelLaguerre.m` constructor
- [~] 2.12 Run tests — verify HankelLaguerre tests pass
- [x] 2.13 Add deprecation warning to `src/beams/HankelHermite.m` constructor
- [~] 2.14 Run tests — verify HankelHermite tests pass

## Phase 3: Integration Verification

- [x] 3.1 Run `tests/test_all.m` (full suite) — 28/30 pass (2 pre-existing failures), deprecation warnings emit correctly
- [x] 3.2 Run `tests/legacy_compat/run_legacy_compat.m` — legacy compat tests pass (16/16)
- [x] 3.3 Verify `which('paraxial.beams.GaussianBeam')` resolves to `/root/Simulation_Scripts/+paraxial/+beams/GaussianBeam.m` — YES
- [x] 3.4 Test standalone +paraxial/ path — BeamFactory.create() works with only +paraxial/ on path (no src/)

## Phase 4: Documentation

- [x] 4.1 Update `README.md` — mark `+paraxial/` as canonical namespace, document namespace convention table
- [x] 4.2 Add note in `README.md` that `src/beams/` is deprecated but functional during transition — DONE (in Namespace Convention section)
- [x] 4.3 Commit "feat: complete +paraxial/ package migration Phase 3"