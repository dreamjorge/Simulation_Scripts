# Consolidation & Modernization Roadmap — Design

**Date:** 2026-04-25
**Status:** Phase A Complete — Phase B In Progress

---

## Integration Result (2026-04-25)

**Master merged with `integration/modernization-v1` (commit `5010c12`)**

| Change | Status |
|--------|--------|
| MATLAB constructor fix (super() first) | ✅ Integrated |
| ~380 test suite | ✅ Integrated |
| OpticalRay / CylindricalRay classes | ✅ Integrated |
| GridUtils asymmetric grids fix | ✅ Integrated |
| HankelLaguerre H1/H2 formula fix | ✅ Integrated |
| CI split portable/legacy | ✅ Integrated |
| Legacy compat layer (16 tests, 16/16 pass) | ✅ Integrated |
| src/ + +paraxial/ structure | ✅ Integrated |
| Documentation alignment | ✅ Integrated |

**Tagged:** `v2.0.0`

---

## Goal

Consolidate valuable changes from `fix/gaussianbeam-super-constructor`, `refactor/utility-classes`, and `chore/legacy-migration-week1` into `master` via selective cherry-pick, and define a clear modernization roadmap toward `+paraxial/` as the canonical namespace.

---

## Part A: Selective Integration

### Source Branches

| Branch | Ahead of master | Role |
|--------|----------------|------|
| `fix/gaussianbeam-super-constructor` | ~254 commits | Primary — contains all valuable changes |
| `refactor/utility-classes` | ~208 commits | Subset of fix/branch — same lineage |
| `chore/legacy-migration-week1` | ~244 commits | Subset of fix/branch — same lineage |

**All three converge on `fix/gaussianbeam-super-constructor`** — it is the parent branch with full history.

### Target Branch

`integration/modernization-v1` created from current `master`.

### Cherry-Pick Groups (in order)

| Group | Key Commit(s) | Description | Verification |
|-------|---------------|-------------|--------------|
| **1. Critical-beam-fixes** | `15a54e9` | Fix MATLAB constructor order (super() unconditional) in all beam subclasses | — |
| **2. Core-tests** | `62acb68`..`aa0f953` | ~380 modular tests with coverage reporting | `octave --no-gui --eval "run('tests/test_all.m')"` |
| **3. Ray-classes** | `4a8b365`, `8b49f8f` | OpticalRay, CylindricalRay classes + tests | — |
| **4. GridUtils-fix** | `bbf5b91` | Asymmetric grids support in static methods | — |
| **5. HankelLaguerre-formula** | `893a397` | Correct H1=LB+i*XLG, H2=LB-i*XLG formula | — |
| **6. CI-split** | `631fec9` | Portable vs legacy compatibility workflow split | CI green |
| **7. Legacy-compat-layer** | From `chore/legacy-migration-week1` | legacy/compat layer, deprecation warnings, smoke runner | legacy_compat suite passes |
| **8. Migration-structure** | `85a19e3`, `1980c29` | src/ directory layout, +paraxial/ Week 7 placeholder | — |
| **9. Docs-alignment** | `96695c9`, `c0be3c1` | Coordinate systems, elegant formulas, canonical examples | — |

### Integration Rules

1. Each group is cherry-picked as a single logical commit where possible
2. After each group: run test suite to verify no regressions
3. If tests fail: fix within the group before proceeding
4. After all groups: full suite run in Octave AND MATLAB before merge to master
5. Merge strategy: `--no-ff merge commit` to preserve logical grouping

### What is NOT integrated (intentionally excluded)

- Anything marked as "WIP" or "in progress" in commit messages
- Partial/incomplete migration phases from Week 1-6 (only Week 7 placeholder included)
- Any changes that break the portable test runner

---

## Part B: Modernization Roadmap

### Vision

The `+paraxial/` package becomes the canonical namespace. Legacy code under `ParaxialBeams/` and `src/` is gradually migrated using the Strangler Fig pattern.

### Phases

```
Post-Integration
│
├── Phase 1: Stabilize
│   ├── All tests green on master after integration
│   ├── Public API contract locked and documented
│   └── No regressions in Octave 11.1.0+ or MATLAB R2020b+
│
├── Phase 2: OO Cleanup (Post-Merge Track from plan.md)
│   ├── Separate model vs computation in Parameters classes
│   ├── Computation layer exists (BeamComputation, HermiteComputation)
│   └── Parameter classes should only hold state + delegate formulas
│
├── Phase 3: +paraxial/ Package Completion
│   ├── Migrate beam classes one-by-one under +paraxial/
│   ├── Maintain backwards compatibility via adapter layer
│   └── Strangler pattern: old code slowly dies as new code grows
│
└── Phase 4: Legacy Deprecation
    ├── Remove Hankele* aliases (already deprecated)
    ├── Clean up legacy/examples/archive
    └── Mark project as deprecated in README
```

### Architectural Principles for +paraxial/

1. **Package structure mirrors responsibility:** `+paraxial/+beams/`, `+paraxial/+parameters/`, `+paraxial/+propagation/+field/`, etc.
2. **No classdef folders** — individual `.m` files per class
3. **IPropagator strategy pattern** remains the propagation interface
4. **Parameters classes** are data/state facades — formula logic stays in computation layer
5. **Backwards compatible** — old `BeamFactory.create('gaussian', ...)` still works via adapters

---

## Decisions

1. **Merge to master:** `git checkout master && git merge --no-ff integration/modernization-v1`
2. **Verification mandatory:** Tests pass in both Octave 11.1.0+ and MATLAB R2020b+ before merge
3. **+paraxial/ placeholder integrated as-is** — skeleton only, no breaking changes
4. **Legacy compat layer integrated complete** — already has tests and smoke runner

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Cherry-pick conflicts on 254 commits | Group by logical blocks; resolve incrementally |
| Tests break during integration | Run after each group; fix before proceeding |
| Legacy compat layer too green | Has smoke runner; will verify before integration |
| MATLAB constructor behavior varies | Critical-beam-fixes group (15a54e9) addresses this explicitly |

---

## Next Step

Invoke **writing-plans skill** to create detailed implementation plan with tasks, checkpoints, and verification criteria.