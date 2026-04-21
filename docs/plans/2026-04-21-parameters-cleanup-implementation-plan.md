# Parameters Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Separate parameter data modeling from beam computation logic to reduce fragility and improve extensibility, without breaking existing public APIs.

**Architecture:** Introduce a new stateless computation layer (`src/computation`) and make parameter classes delegate numerical formulas to it. Keep all current constructors/methods/properties stable, including legacy shims, then verify behavior with focused formula tests + existing portable/legacy suites.

**Tech Stack:** MATLAB/Octave (`.m` classdef + script tests), existing `tests/test_all.m` + `tests/portable_runner.m` + legacy compatibility runner.

---

### Task 1: Create stateless BeamComputation foundation

**Files:**
- Create: `src/computation/BeamComputation.m`
- Create: `tests/modern/test_BeamComputation.m`

**Step 1: Write the failing test file for core formulas**

Create `tests/modern/test_BeamComputation.m` with checks for:
- `rayleighDistance(w0, lambda)` equals `pi*w0^2/lambda`
- `waveNumber(lambda)` equals `2*pi/lambda`
- `waist(w0, 0, lambda, zr) == w0`
- `waist(w0, zr, lambda, zr) == w0*sqrt(2)`
- `radiusOfCurvature(0, zr)` is `Inf`
- `complexBeamParameter(z, zr, k) == 1i*k/(2*(z+1i*zr))`

**Step 2: Run test to verify it fails**

Run:  
`& "C:\Users\uidn7961\AppData\Local\Programs\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe" --no-gui --eval "run('tests/modern/test_BeamComputation.m')"`

Expected: FAIL (class/function `BeamComputation` not found).

**Step 3: Implement BeamComputation minimal API**

Create `src/computation/BeamComputation.m` with static methods:
- `rayleighDistance`
- `waveNumber`
- `waist`
- `gouyPhase`
- `radiusOfCurvature`
- `complexBeamParameter`

Use vectorized math and Octave-safe operations.

**Step 4: Run test to verify it passes**

Run same command from Step 2.

Expected: PASS for all core formula checks.

**Step 5: Commit**

```bash
git add src/computation/BeamComputation.m tests/modern/test_BeamComputation.m
git commit -m "feat(parameters): add stateless beam computation utilities"
```

---

### Task 2: Wire GaussianParameters to delegation layer

**Files:**
- Modify: `src/parameters/GaussianParameters.m`
- Test: `tests/modern/test_GaussianParameters.m`

**Step 1: Add/adjust failing assertions for delegation-equivalent behavior**

In `tests/modern/test_GaussianParameters.m`, add checks that compare:
- `params.waist(z2)` against direct formula (`BeamComputation.waist(...)`)
- `params.gouyPhase(z2)` against `BeamComputation.gouyPhase(...)`
- `params.radius(z2)` against `BeamComputation.radiusOfCurvature(...)`

Do not remove legacy snapshot checks (`Waist`, `GouyPhase`, `Radius`).

**Step 2: Run GaussianParameters test to verify current status**

Run:  
`& "C:\Users\uidn7961\AppData\Local\Programs\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe" --no-gui --eval "run('tests/modern/test_GaussianParameters.m')"`

Expected: If new assertions reference BeamComputation path before setpaths updates, may FAIL.

**Step 3: Refactor GaussianParameters to delegate formulas**

Update constructor and dynamic methods to call `BeamComputation` for formula computation while keeping:
- constructor signature unchanged
- property names unchanged
- dynamic method names unchanged
- static legacy methods unchanged behavior

**Step 4: Run GaussianParameters test again**

Run same command from Step 2.

Expected: PASS with no behavior regression.

**Step 5: Commit**

```bash
git add src/parameters/GaussianParameters.m tests/modern/test_GaussianParameters.m
git commit -m "refactor(parameters): delegate gaussian formulas to BeamComputation"
```

---

### Task 3: Delegate Hermite/Laguerre parameter families

**Files:**
- Modify: `src/parameters/HermiteParameters.m`
- Modify: `src/parameters/LaguerreParameters.m`
- Modify: `src/parameters/ElegantHermiteParameters.m`
- Modify: `src/parameters/ElegantLaguerreParameters.m`
- Test: `tests/modern/test_HermiteParameters.m`
- Test: `tests/modern/test_LaguerreParameters.m`
- Test: `tests/modern/test_ElegantHermiteParameters.m`
- Test: `tests/modern/test_ElegantLaguerreParameters.m`

**Step 1: Add failing equivalence checks per family**

For each test file, add assertions validating:
- modal Gouy (`phiPhase`) formula equivalence
- waist scaling equivalence (Hermite/Laguerre)
- elegant alpha equivalence (`alphaAtZ` and `alpha`)

**Step 2: Run each parameter-family test to capture failures**

Run each test independently with octave-cli and confirm failure mode.

**Step 3: Refactor each parameter class to delegate formulas**

Use `BeamComputation` methods instead of duplicating formulas inline.
Preserve all existing public methods/properties and default values.

**Step 4: Re-run each parameter-family test**

Expected: all updated tests PASS with old behavior preserved.

**Step 5: Commit**

```bash
git add src/parameters/HermiteParameters.m src/parameters/LaguerreParameters.m src/parameters/ElegantHermiteParameters.m src/parameters/ElegantLaguerreParameters.m tests/modern/test_HermiteParameters.m tests/modern/test_LaguerreParameters.m tests/modern/test_ElegantHermiteParameters.m tests/modern/test_ElegantLaguerreParameters.m
git commit -m "refactor(parameters): delegate higher-order families to computation layer"
```

---

### Task 4: Move legacy Hermite helper into HermiteComputation with shim

**Files:**
- Create: `src/computation/HermiteComputation.m`
- Modify: `src/parameters/HermiteParameters.m`
- Test: `tests/legacy_compat/test_HankelCompatibility.m`

**Step 1: Add failing test coverage for helper compatibility path**

Add a test assertion that calls both:
- `HermiteParameters.getHermiteSolutions(nu, x)`
- `HermiteComputation.hermiteSolutions(nu, x)`

and verifies identical outputs within tolerance.

**Step 2: Run compatibility test to verify failure**

Run specific legacy compatibility script (or targeted script section) and confirm missing class/method fails.

**Step 3: Implement HermiteComputation and compatibility shim**

- Move algorithm body into `HermiteComputation.hermiteSolutions`
- Keep `HermiteParameters.getHermiteSolutions` delegating to new class
- Add comment marking old entrypoint as compatibility shim

**Step 4: Re-run compatibility test**

Expected: PASS, with old entrypoint still working.

**Step 5: Commit**

```bash
git add src/computation/HermiteComputation.m src/parameters/HermiteParameters.m tests/legacy_compat/test_HankelCompatibility.m
git commit -m "refactor(parameters): extract Hermite helper with legacy shim"
```

---

### Task 5: Ensure path wiring includes new computation folder

**Files:**
- Modify: `setpaths.m`
- Modify: `tests/portable_runner.m`

**Step 1: Add failing check in portable context**

Ensure runner can resolve `BeamComputation` without manual addpath from CWD.

**Step 2: Run portable suite to verify failure (if path missing)**

Run:  
`& "C:\Users\uidn7961\AppData\Local\Programs\GNU Octave\Octave-11.1.0\mingw64\bin\octave-cli.exe" --no-gui --eval "run('tests/test_all.m')"`

Expected: failure if `src/computation` not on path.

**Step 3: Update path setup**

- Add `src/computation` to `setpaths.m`
- Add `src/computation` to `tests/portable_runner.m`

**Step 4: Re-run portable suite**

Expected: suite completes successfully.

**Step 5: Commit**

```bash
git add setpaths.m tests/portable_runner.m
git commit -m "chore(paths): include computation layer in setup and portable runner"
```

---

### Task 6: Update architecture docs and verify full regression

**Files:**
- Modify: `docs/ARCHITECTURE.md`
- Modify: `README.md` (only if needed for structure mention)

**Step 1: Add failing doc checklist (manual)**

Checklist:
- mentions `src/computation` layer
- shows dependency direction `Parameters -> BeamComputation`
- clarifies legacy shim policy

**Step 2: Update architecture doc minimally**

Document:
- new computation layer purpose
- unchanged public APIs
- migration rationale (fragility/extensibility)

**Step 3: Run complete regression suites**

Run in order:

1. `tests/modern/test_BeamComputation.m`
2. `tests/test_all.m`
3. `tests/legacy_compat/run_legacy_compat.m`

Expected: all green.

**Step 4: Final API compatibility smoke checks**

Run canonical examples (non-build execution) as quick smoke:
- `examples/canonical/MainGauss_refactored.m`
- `examples/canonical/MainMultiMode.m`
- `examples/canonical/ExampleRayTracing.m`

Expected: no signature/runtime regressions related to parameter access.

**Step 5: Commit**

```bash
git add docs/ARCHITECTURE.md README.md
git commit -m "docs(architecture): document computation layer and delegation model"
```

---

### Task 7: Final readiness gate for merge

**Files:**
- Modify: `plan.md` (optional status note)

**Step 1: Verify no public API drift**

Confirm unchanged signatures for:
- `ParaxialBeam` contract (`opticalField`, `getParameters`, `beamName`)
- all parameter constructors
- `BeamFactory.create`

**Step 2: Verify constraints still hold**

- Octave/MATLAB compatibility maintained
- no external dependencies added
- tests green

**Step 3: Produce short verification summary**

Record:
- changed files
- test commands executed
- pass/fail status
- known deferred items (if any)

**Step 4: Commit final status note (if `plan.md` updated)**

```bash
git add plan.md
git commit -m "docs: record parameters cleanup verification status"
```

---

## Definition of Done

- `src/computation/BeamComputation.m` exists and is covered by tests.
- parameter classes delegate formulas without breaking public behavior.
- Hermite legacy helper extracted with compatibility shim preserved.
- `setpaths.m` and `portable_runner.m` include `src/computation`.
- portable + legacy tests pass.
- architecture docs reflect the new dependency direction.
