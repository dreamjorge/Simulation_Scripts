# Proposal: Post-Merge Cleanup — After v2.0.0 Integration

## Intent

Clean up remaining technical debt after v2.0.0 integration: register missing test, prune stale branches, formalize legacy policy, and advance the Strangler Fig migration from `src/` to `+paraxial/`.

## Scope

### In Scope
- Register `test_Wavefront.m` in `test_all.m` (1-line fix, task 8.3 of wavefront-class)
- Delete stale branches: `chore/legacy-migration-week1`, `fix/gaussianbeam-super-constructor`, `refactor/utility-classes`
- Create `LEGACY_POLICY.md` defining what to do with 33 examples in `examples/legacy/`
- Emit deprecation warnings in `src/` classes (already done for beams; extend to `src/propagation/`)
- Archive completed SDD changes

### Out of Scope
- Running performance benchmarks (blocked: no MATLAB license / Octave runtime in CI)
- Deep OO redesign
- Package migration completion (requires coordinated example updates)

## Capabilities

### New Capabilities
- `legacy-policy`: Documented policy for managing 33 legacy examples in `examples/legacy/`
- `deprecation-propagation`: Extend deprecation warnings from `src/beams/` to `src/propagation/` classes

### Modified Capabilities
- None

## Approach

1. **Test registration**: Add `test_Wavefront.m` to `test_all.m` test list
2. **Branch cleanup**: Delete stale local and remote branches via `git branch -d` and `git push origin --delete`
3. **Legacy policy**: Write `examples/legacy/LEGACY_POLICY.md` if not exists, else verify it matches current state
4. **Deprecation extension**: Add deprecation warnings to `src/propagation/field/` and `src/propagation/rays/` classes
5. **Archive completed SDDs**: Update `state.yaml` for completed changes to `completed` status

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `tests/test_all.m` | Modified | Register test_Wavefront.m |
| `examples/legacy/LEGACY_POLICY.md` | Modified | Clarify retention/removal policy |
| `src/propagation/field/` | Modified | Add deprecation warnings |
| `src/propagation/rays/` | Modified | Add deprecation warnings |
| Git branches | Removed | Delete 3 stale branches |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Branch deletion causes merge conflicts later | Low | All stale branches are 200+ commits behind master |
| Legacy policy contradicts actual usage | Low | Verify against actual usage signals before writing |

## Rollback Plan

```bash
# Revert test registration:
git checkout HEAD~1 -- tests/test_all.m

# Re-add branches:
git branch chore/legacy-migration-week1 <sha>
git branch fix/gaussianbeam-super-constructor <sha>
git branch refactor/utility-classes <sha>

# Revert deprecation warnings:
git checkout HEAD~1 -- src/propagation/field/*.m src/propagation/rays/*.m
```

## Dependencies

- GitHub Actions CI passing for Octave and MATLAB

## Success Criteria

- [ ] `test_Wavefront.m` runs in `portable_runner` suite
- [ ] 3 stale branches deleted locally and remotely
- [ ] `examples/legacy/LEGACY_POLICY.md` exists and reflects current classification (archive/generators/research)
- [ ] `src/propagation/field/IPropagator.m`, `FFTPropagator.m`, `AnalyticPropagator.m` emit deprecation warnings
- [ ] `src/propagation/rays/` classes emit deprecation warnings
- [ ] Completed SDD changes have `state.yaml` updated to `completed`
