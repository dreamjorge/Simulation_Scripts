# Proposal: +paraxial/ Package Migration (Phase 3)

## Intent

Unify the duplicate codebase structure by making `+paraxial/` the canonical namespace for beam classes. Currently `src/` and `+paraxial/` are parallel copies of nearly identical code. This creates maintenance confusion and violates single-source-of-truth. The Strangler Fig pattern will gradually move all canonical implementation to `+paraxial/`, deprecating `src/` paths.

## Scope

### In Scope
- Migrate `src/beams/*.m` → `+paraxial/+beams/*.m` (beam classes one-by-one)
- Migrate `src/parameters/*.m` → `+paraxial/+parameters/*.m` (parameter classes)
- Update `ParaxialBeams/BeamFactory.m` to prefer `+paraxial/` classes
- Add deprecation warnings to `src/` classes after migration
- Verify full test suite passes after each class migration
- Ensure `+paraxial/` classes work standalone (no cross-namespace dependencies)

### Out of Scope
- Migrating `src/propagation/` — already has partial `+paraxial/+propagation/+rays/`
- Migrating `src/computation/` — stateless utility, used by both namespaces
- Migrating `src/visualization/` — separate concern, future work
- Breaking existing `src/` API contracts during transition

## Capabilities

### New Capabilities
- `paraxial-beam-factory`: BeamFactory creates instances from `+paraxial/` classes instead of `src/`

### Modified Capabilities
- None — this is a pure refactor (code move + namespace consolidation)

## Approach

**Strangler Fig pattern** — `+paraxial/` becomes canonical, `src/` becomes adapter layer.

```
Sequence per beam class:
  1. Verify +paraxial/ version exists and tests pass
  2. Update BeamFactory to instantiate from +paraxial/
  3. Add @deprecation warning to src/ version
  4. Run full test suite — must pass
  5. Commit

Order: GaussianBeam → HermiteBeam → LaguerreBeam → ElegantHermiteBeam →
       ElegantLaguerreBeam → HankelLaguerre → HankelHermite
```

After all beams migrated: deprecate `src/beams/` directory with notice in README.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `+paraxial/+beams/*.m` | Modified | New canonical location |
| `+paraxial/+parameters/*.m` | Modified | New canonical location |
| `src/beams/*.m` | Deprecated | Becomes thin adapter/wrapper |
| `src/parameters/*.m` | Deprecated | Becomes thin adapter/wrapper |
| `ParaxialBeams/BeamFactory.m` | Modified | Routes to +paraxial/ classes |
| `tests/test_all.m` | Verified | Ensures no regressions |
| `README.md` | Modified | Documents new canonical paths |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Breaking existing code that imports `src/beams/` directly | Medium | Keep src/ as functional adapter during transition |
| `+paraxial/` classes have hidden `src/` dependencies | Low | Test each class in isolation before updating BeamFactory |
| Octave/MATLAB package path differences | Low | Use `which()` to verify class resolution before/after |

## Rollback Plan

If migration breaks tests: revert BeamFactory change, remove deprecation warnings from src/ classes, and revert README.md changes. The src/ classes remain functional until fully cut over.

## Dependencies

- `@deprecation` helper function (MATLAB-only) — for Octave, use warning() instead

## Success Criteria

- [ ] BeamFactory creates beams from `+paraxial/+beams/` classes
- [ ] All `tests/test_all.m` tests pass (Octave + MATLAB)
- [ ] All `tests/legacy_compat/run_legacy_compat.m` pass
- [ ] `src/beams/` classes have deprecation warnings
- [ ] README.md documents `+paraxial/` as canonical
- [ ] No cross-namespace dependencies in `+paraxial/` classes (verified via `which()`)