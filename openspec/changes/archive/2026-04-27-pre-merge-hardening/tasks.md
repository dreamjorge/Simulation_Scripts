# Tasks: Pre-Merge Hardening Documentation

## Phase 1: Read Existing Code and Verify State

- [ ] 1.1 Read `README.md` current content and list discrepancies with actual repo
- [ ] 1.2 List all 27 files in `ParaxialBeams/*.m` for reference
- [ ] 1.3 List all files in `examples/*.m` for classification
- [ ] 1.4 Compare `exec/pre-merge-hardening` vs `integration/pre-master` to identify redundant work
- [ ] 1.5 Read current `plan.md` pre-merge checklist items

## Phase 2: Create Architecture Documentation

- [ ] 2.1 Create `docs/ARCHITECTURE.md` with class hierarchy diagram
- [ ] 2.2 Document Strategy Pattern (IPropagator interface and 3 implementations)
- [ ] 2.3 Document Factory Pattern (BeamFactory.create)
- [ ] 2.4 Document Beam API contract (opticalField, getParameters, beamName)
- [ ] 2.5 Document data flow: Grid → Beam → Propagator → Result

## Phase 3: Update README.md

- [ ] 3.1 Rewrite structure section to match actual 27 .m files in ParaxialBeams/
- [ ] 3.2 Update class descriptions (PhysicalConstants, GridUtils, FFTUtils, etc.)
- [ ] 3.3 Add canonical examples section listing MainGauss_refactored.m, MainMultiMode.m, ExampleRayTracing.m
- [ ] 3.4 Update MATLAB/Octave compatibility section
- [ ] 3.5 Add Beam API Contract section documenting the 3-method interface
- [ ] 3.6 Verify all file paths in README match actual repo structure

## Phase 4: Classify Examples

- [ ] 4.1 Add `%% canonical` header to examples/MainGauss_refactored.m
- [ ] 4.2 Add `%% canonical` header to examples/MainMultiMode.m
- [ ] 4.3 Add `%% canonical` header to ExampleRayTracing.m
- [ ] 4.4 Add `%% legacy` header to remaining examples/ scripts
- [ ] 4.5 Run canonical examples to verify they execute without error

## Phase 5: Complete plan.md Checklist

- [ ] 5.1 Mark pre-merge scope items as complete in plan.md
- [ ] 5.2 Mark API audit items as complete if verified
- [ ] 5.3 Mark canonical examples classification as complete
- [ ] 5.4 Mark test coverage gates as verified
- [ ] 5.5 Add section documenting deferred post-merge work

## Phase 6: Branch Hygiene

- [ ] 6.1 Archive or merge `exec/pre-merge-hardening` branch (decide based on Phase 1.4)
- [ ] 6.2 Verify no conflicts between integration/pre-master and master
- [ ] 6.3 Draft merge commit message following conventional commits
- [ ] 6.4 Create CHANGELOG.md or add to existing documenting changes

## Phase 7: Final Verification

- [ ] 7.1 Run `octave tests/test_all.m` to confirm suite passes
- [ ] 7.2 Verify README.md renders correctly with `cat README.md`
- [ ] 7.3 Confirm docs/ARCHITECTURE.md exists and contains diagrams
- [ ] 7.4 Review final state with `git status` and `git diff master`
