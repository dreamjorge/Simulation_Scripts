# Consolidation & Modernization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate select changes from `fix/gaussianbeam-super-constructor` / `refactor/utility-classes` / `chore/legacy-migration-week1` into `master` via themed cherry-picks, then establish a clear modernization roadmap toward `+paraxial/` as canonical namespace.

**Architecture:** Selective cherry-pick grouped by theme (beam fixes, tests, ray classes, CI, legacy compat, migration structure). Each group verified with test suite before proceeding. Post-integration modernization phased toward `+paraxial/` package.

**Tech Stack:** MATLAB, Octave, `.m` scripts, Git, GitHub Actions, `tests/test_all.m`, `tests/portable_runner.m`.

---

## PHASE A: INTEGRATION

---

### Task 1: Create integration branch

**Files:**
- Reference: `master`

**Step 1: Verify master is clean and up-to-date**

Run: `git status`
Expected: `On branch master. Your branch is up to date with 'origin/master'. nothing to commit, working tree clean.`

**Step 2: Create integration branch from master**

Run: `git checkout -b integration/modernization-v1`
Expected: Switched to a new branch 'integration/modernization-v1'

**Step 3: Verify branch is clean**

Run: `git status`
Expected: `On branch integration/modernization-v1. nothing to commit, working tree clean.`

**Step 4: Push branch to remote**

Run: `git push -u origin integration/modernization-v1`
Expected: Branch created and tracked.

---

### Task 2: Cherry-pick Group 1 — Critical Beam Constructor Fixes

**Files:**
- Modify: `src/beams/GaussianBeam.m`
- Modify: `src/beams/HermiteBeam.m`
- Modify: `src/beams/LaguerreBeam.m`
- Modify: `src/beams/ElegantHermiteBeam.m`
- Modify: `src/beams/ElegantLaguerreBeam.m`
- Modify: `src/beams/HankelLaguerre.m`
- Reference: `src/beams/ParaxialBeam.m`

**Step 1: Cherry-pick the constructor fix commit**

Run: `git cherry-pick 15a54e9`
Expected: `[integration/modernization-v1 1234567] fix(beams): call ParaxialBeam constructor first in all beam subclasses`

**Step 2: Verify cherry-pick succeeded**

Run: `git log --oneline -3`
Expected: Top commit is 15a54e9 with the constructor fix message.

**Step 3: Run tests to verify no regressions**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass. Exit code 0.

**Step 4: Commit (already done by cherry-pick)**

No action needed. Cherry-pick creates its own commit.

---

### Task 3: Cherry-pick Group 2 — Core Test Suite Expansion

**Files:**
- Reference: `tests/modern/test_*.m` (~20 test files)
- Reference: `tests/test_all.m`
- Reference: `tests/runAllTests.m`

**Step 1: Identify test commit range**

Run: `git log fix/gaussianbeam-super-constructor --oneline | grep -E "test: add|test: expand|test: modularize" | head -20`
Expected: List of test commits from `62acb68` through `aa0f953`.

**Step 2: Cherry-pick all test expansion commits as a range**

Run: `git cherry-pick 62acb68..aa0f953`
Expected: Multiple commits applied. Check with `git log --oneline` to see the applied commits.

**Step 3: Run full test suite**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All ~380 tests pass. Exit code 0.

**Step 4: If any test fails, identify the culprit commit and fix it**

Run: `git bisect start` (if needed) to find which commit broke tests.
Fix the issue in a new commit on the integration branch.
Run tests again to verify.

**Step 5: Push changes**

Run: `git push`
Expected: Remote updated with test commits.

---

### Task 4: Cherry-pick Group 3 — Ray Classes (OpticalRay, CylindricalRay)

**Files:**
- Create: `src/propagation/rays/OpticalRay.m`
- Create: `src/propagation/rays/CylindricalRay.m`
- Create: `tests/modern/test_CylindricalRay.m`
- Create: `tests/modern/test_OpticalRay.m`
- Modify: `tests/test_all.m`

**Step 1: Cherry-pick OpticalRay commit**

Run: `git cherry-pick 4a8b365`
Expected: `[integration/modernization-v1...] feat: add OpticalRay class for Cartesian ray tracing`

**Step 2: Cherry-pick CylindricalRay + tests commit**

Run: `git cherry-pick 8b49f8f`
Expected: `[integration/modernization-v1...] test: add CylindricalRay tests and register in test suite`

**Step 3: Verify files exist**

Run: `ls -la src/propagation/rays/OpticalRay.m src/propagation/rays/CylindricalRay.m`
Expected: Both files exist.

**Step 4: Run tests**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass including new ray tests.

**Step 5: Push**

Run: `git push`
Expected: Remote updated.

---

### Task 5: Cherry-pick Group 4 — GridUtils Asymmetric Fix

**Files:**
- Modify: `ParaxialBeams/GridUtils.m`

**Step 1: Cherry-pick GridUtils fix**

Run: `git cherry-pick bbf5b91`
Expected: `[integration/modernization-v1...] fix(GridUtils): support asymmetric grids in static methods`

**Step 2: Run tests**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass.

**Step 3: Push**

Run: `git push`

---

### Task 6: Cherry-pick Group 5 — HankelLaguerre Formula Fix

**Files:**
- Modify: `src/beams/HankelLaguerre.m`

**Step 1: Cherry-pick HankelLaguerre formula commit**

Run: `git cherry-pick 893a397`
Expected: `[integration/modernization-v1...] feat(ParaxialBeams): implement HankelLaguerre with H1=LB+i*XLG, H2=LB-i*XLG formula`

**Step 2: Verify HankelLaguerre tests still pass**

Run: `octave --no-gui --eval "run('tests/modern/test_HankelLaguerre.m')"`
Expected: All HankelLaguerre tests pass.

**Step 3: Run full suite**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass.

**Step 4: Push**

Run: `git push`

---

### Task 7: Cherry-pick Group 6 — CI Split Portable/Legacy

**Files:**
- Modify: `.github/workflows/octave.yml`
- Modify: `.github/workflows/matlab.yml`

**Step 1: Cherry-pick CI split commit**

Run: `git cherry-pick 631fec9`
Expected: `[integration/modernization-v1...] ci(workflows): split portable and legacy compatibility suites`

**Step 2: Verify workflows have both portable and legacy jobs**

Run: `cat .github/workflows/octave.yml | grep -E "portable|legacy"`
Expected: Workflow file contains both portable and legacy compatibility job definitions.

**Step 3: Push**

Run: `git push`

---

### Task 8: Cherry-pick Group 7 — Legacy Compat Layer

**Files:**
- Create: `legacy/compat/` directory and files
- Create: `tests/legacy_compat/*.m`
- Modify: `tests/run_legacy_compat.m`

**Step 1: Cherry-pick legacy compat commits from chore/legacy-migration-week1**

Run: `git cherry-pick 3c24073 175eec9 2735981 7df49ec`
Expected: Legacy compat layer commits applied.

**Step 2: Verify legacy compat tests exist**

Run: `ls tests/legacy_compat/`
Expected: `test_HankelAliasEdgeCases.m`, `test_HankelAliasStaticDelegation.m`, `test_HankelCompatibility.m`, `test_LegacyBeamConstructors.m`.

**Step 3: Run legacy compat suite**

Run: `octave --no-gui --eval "run('tests/legacy_compat/run_legacy_compat.m')"`
Expected: All legacy compat tests pass.

**Step 4: Push**

Run: `git push`

---

### Task 9: Cherry-pick Group 8 — Migration Structure (src/ layout + +paraxial/ placeholder)

**Files:**
- Create: `src/beams/`, `src/parameters/`, `src/propagation/`, `src/computation/`, `src/visualization/`
- Modify: `+paraxial/` (Week 7 placeholder)

**Step 1: Cherry-pick src/ directory layout commit**

Run: `git cherry-pick 85a19e3`
Expected: `[integration/modernization-v1...] refactor(structure): implement src/ directory layout per migration plan`

**Step 2: Cherry-pick +paraxial/ placeholder commit**

Run: `git cherry-pick 1980c29`
Expected: `[integration/modernization-v1...] refactor(migration): implement Week 7 +paraxial/ package placeholder`

**Step 3: Verify +paraxial/ placeholder exists**

Run: `ls +paraxial/`
Expected: Package structure exists (may be minimal/placeholder).

**Step 4: Run tests**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass.

**Step 5: Push**

Run: `git push`

---

### Task 10: Cherry-pick Group 9 — Documentation Alignment

**Files:**
- Modify: `docs/ARCHITECTURE.md`
- Modify: `README.md`

**Step 1: Cherry-pick coordinate systems docs**

Run: `git cherry-pick 96695c9`
Expected: `[integration/modernization-v1...] docs: document coordinate systems and elegant beam formulas`

**Step 2: Cherry-pick canonical examples docs**

Run: `git cherry-pick c0be3c1`
Expected: `[integration/modernization-v1...] docs: update documentation for new src/ structure and Week 3 progress`

**Step 3: Verify docs are consistent**

Run: `grep -l "canonical" README.md docs/ARCHITECTURE.md`
Expected: Both files reference canonical examples correctly.

**Step 4: Push**

Run: `git push`

---

### Task 11: Final Verification — Full Test Suite in Octave AND MATLAB

**Files:**
- Reference: All source files
- Reference: All test files

**Step 1: Run full test suite in Octave**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass. Exit code 0.

**Step 2: Run legacy compat suite in Octave**

Run: `octave --no-gui --eval "run('tests/legacy_compat/run_legacy_compat.m')"`
Expected: All legacy compat tests pass.

**Step 3: Verify CI workflows are configured correctly**

Run: `cat .github/workflows/octave.yml`
Expected: Workflow contains portable AND legacy jobs, triggers on push to integration/*.

**Step 4: Push final state**

Run: `git push`
Expected: All commits pushed.

---

### Task 12: Merge to Master

**Files:**
- Modify: `master`

**Step 1: Verify all tasks complete**

Run: `git log integration/modernization-v1 --oneline | wc -l`
Expected: All cherry-picked commits present (at least 15+ commits).

**Step 2: Run final full test suite**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass before merge.

**Step 3: Checkout master**

Run: `git checkout master`
Expected: Switched to branch 'master'.

**Step 4: Merge with --no-ff to preserve logical grouping**

Run: `git merge --no-ff integration/modernization-v1`
Expected: Merge commit created successfully.

**Step 5: Run tests on master after merge**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: All tests pass on master.

**Step 6: Push master**

Run: `git push`
Expected: Master updated on remote.

**Step 7: Tag the release**

Run: `git tag -a v2.0.0 -m "Consolidation: beam fixes, ~380 tests, ray classes, CI split, legacy compat, +paraxial/ placeholder"`
Expected: Tag created.

**Step 8: Push tags**

Run: `git push origin v2.0.0`
Expected: Tag pushed to remote.

---

## PHASE B: MODERNIZATION ROADMAP

---

### Task 13: Document Post-Integration Modernization Phases

**Files:**
- Modify: `docs/plans/2026-04-25-consolidation-modernization-design.md`

**Step 1: Add Phase 1 details to design doc**

Add to design doc:

```markdown
### Phase 1: Stabilize Integration
- Lock public API contract (already documented in README.md)
- Verify no regressions in Octave 11.1.0+ and MATLAB R2020b+
- All critical paths tested

### Phase 2: OO Cleanup
- Parameter classes should only hold state (w0, lambda, z) + delegate to computation layer
- BeamComputation and HermiteComputation already exist as formula layer
- Clean up any inline formula logic remaining in parameter classes

### Phase 3: +paraxial/ Package Completion
- Migrate one beam class at a time under +paraxial/
- Keep adapters so BeamFactory.create() still works
- Test after each migration

### Phase 4: Legacy Deprecation
- Remove Hankele* aliases (already deprecated)
- Clean legacy/examples/archive
- Mark deprecated in README
```

**Step 2: Commit Phase B roadmap**

Run: `git add docs/plans/2026-04-25-consolidation-modernization-design.md && git commit -m "docs: add modernization phases to design doc"`

---

## COMPLETION

**Plan complete and saved to `docs/plans/2026-04-25-consolidation-modernization-implementation-plan.md`.**

Two execution options:

**1. Subagent-Driven (this session)** — I dispatch fresh subagent per task, review between tasks, fast iteration

**2. Parallel Session (separate)** — Open new session with executing-plans, batch execution with checkpoints

**Which approach?**