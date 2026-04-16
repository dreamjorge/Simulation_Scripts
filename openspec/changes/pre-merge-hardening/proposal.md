# Proposal: Pre-Merge Hardening — `integration/pre-master` → `master`

## Intent

Consolidate and document the current state of `integration/pre-master` for a clean merge to `master`. The branch has 71 commits of massive refactoring (Strategy/Factory patterns, ~380 tests, CI/CD for MATLAB/Octave) but lacks documentation that reflects the reality of the code. Without this work, the merge perpetuates the misalignment between docs and code.

## Scope

### In Scope
- Update `README.md` with the actual repo structure (27 .m files in ParaxialBeams)
- Create `docs/ARCHITECTURE.md` with class diagram and applied patterns
- Document the public beam API (contract: `opticalField`, `getParameters`, `beamName`)
- Classify examples as `canonical` vs `legacy` in `README.md`
- Complete the merge readiness checklist in `plan.md`
- Add CHANGELOG.md or conventional commits for the merge
- Unify redundant branches (`exec/pre-merge-hardening` vs `integration/pre-master`)

### Out of Scope
- Package migration to `+paraxial/` namespaces (post-merge)
- Deep OO redesign of beams/propagators (post-merge)
- Full rewrite of historical examples (post-merge)
- Rescue of complete legacy branches
- Major structural cleanup of third-party addons

## Capabilities

### New Capabilities
- `pre-merge-docs`: Architecture and API documentation to facilitate merge to master

### Modified Capabilities
- `beam-api`: The existing ParaxialBeam contract is formally documented but does not change

## Approach

1. **API Audit**: Read the 6 beam class files and verify they comply with the `opticalField(X,Y,z)`, `getParameters(z)`, `beamName()` contract
2. **Canonical Documentation**: Create `docs/ARCHITECTURE.md` with class structure, patterns, and data flow
3. **README.md rewrite**: Replace old structure (`@Folder/`) with actual reality (`.m` files)
4. **Example classification**: Mark `examples/MainGauss_refactored.m`, `examples/MainMultiMode.m`, `ExampleRayTracing.m` as canonical
5. **Branch unification**: Compare `exec/pre-merge-hardening` with `integration/pre-master` and consolidate or archive

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `README.md` | Modified | Actual structure, canonical examples, MATLAB/Octave compatibility |
| `plan.md` | Modified | Complete readiness checklist for pre-merge |
| `docs/ARCHITECTURE.md` | New | Class diagram, Strategy/Factory patterns, data flow |
| `examples/` | Modified | Classify scripts as canonical vs legacy |
| `openspec/` | New | SDD structure for tracking |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Conflicts between `exec/pre-merge-hardening` and `integration/pre-master` | Medium | Diff both against master and decide which to absorb |
| Outdated `README.md` confuses users post-merge | High | Cross-verify with `ls ParaxialBeams/*.m` before commit |
| Broken tests in MATLAB (only verified in Octave) | Low | CI already has MATLAB workflow |

## Rollback Plan

```bash
# If the merge has issues:
git revert <merge-commit>
git checkout master
# Restore docs from the commit prior to the merge
```

## Dependencies

- GitHub Actions CI passing for Octave and MATLAB (already configured)
- 71 refactor commits in `integration/pre-master` verified

## Success Criteria

- [ ] `README.md` reflects actual structure: `ls ParaxialBeams/*.m | wc -l` = 27 files documented
- [ ] `docs/ARCHITECTURE.md` exists with architecture diagram
- [ ] `plan.md` has all pre-merge checklist items marked as complete
- [ ] Canonical examples identified and documented
- [ ] `exec/pre-merge-hardening` archived or merged
- [ ] Merge commit follows conventional commits: `merge(integration): stabilize beam API...`
