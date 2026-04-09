# Integration Pre-Merge Hardening Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** dejar `integration/pre-master` listo para mergearse a `master/main` con alcance congelado, API pública auditada, ejemplos canónicos definidos y evidencia real de compatibilidad MATLAB/Octave.

**Architecture:** este trabajo es de hardening, no de rediseño. La implementación se apoya en tres superficies reales del repo: código fuente en `ParaxialBeams/`, pruebas en `tests/` y narrativa pública en `README.md`. Cada cambio tiene que reducir ambiguedad sin reabrir la arquitectura OO.

**Tech Stack:** MATLAB, Octave, scripts `.m`, test runner `tests/test_all.m`, runner portátil `tests/portable_runner.m`, Git/GitHub Actions.

## Merge Scope Checklist

- [ ] `ParaxialBeams/*.m` modernos auditados
- [ ] `tests/` vigentes y ejecutables en Octave/MATLAB
- [ ] compatibilidad y portabilidad preservadas
- [ ] ejemplos canónicos identificados
- [ ] narrativa pública alineada con el estado real

## Explicitly Out of Scope

- package migration a `+paraxial/...`
- rediseño OO profundo de beams/propagators
- reescritura total de ejemplos históricos
- rescate de ramas legacy completas
- limpieza estructural grande de addons terceros

---

### Task 1: Freeze Merge Scope

**Files:**
- Modify: `plan.md`
- Reference: `README.md`

**Step 1: Document the in-scope deliverables**

Agregar una checklist de merge scope en `plan.md` con estos items exactos:

```md
- [ ] `ParaxialBeams/*.m` modernos auditados
- [ ] `tests/` vigentes y ejecutables en Octave/MATLAB
- [ ] compatibilidad y portabilidad preservadas
- [ ] ejemplos canónicos identificados
- [ ] narrativa pública alineada con el estado real
```

**Step 2: Document the out-of-scope work**

Agregar una sección `## Explicitly Out of Scope` en `plan.md` con estos items:

```md
- package migration a `+paraxial/...`
- rediseño OO profundo de beams/propagators
- reescritura total de ejemplos históricos
- rescate de ramas legacy completas
- limpieza estructural grande de addons terceros
```

**Step 3: Verify the plan still describes a stabilization merge**

Run: revisar `plan.md` y confirmar que no aparece trabajo de rediseño como prerequisito del merge.
Expected: el documento separa explícitamente hardening pre-merge vs. rediseño post-merge.

**Step 4: Commit**

```bash
git add plan.md
git commit -m "docs: freeze pre-merge scope"
```

### Task 2: Audit Public Beam API

| Class | Primary field method | Accepted coordinates | z semantics | Uses Parameters state | Status |
| --- | --- | --- | --- | --- | --- |
| `ParaxialBeam` | `opticalField(obj, X, Y, z)` | cartesian `X`, `Y` | abstract contract requires dynamic `z` argument | no | aligned |
| `GaussianBeam` | `opticalField(obj, X, Y, z)` | cartesian `X`, `Y` converted to radial distance internally | dynamic `z` used directly to compute `w`, `R`, and `psi` | yes, for `InitialWaist` and `Lambda` only | aligned |
| `HermiteBeam` | `opticalField(obj, X, Y, z)` -> `computeComplexField(X, Y, obj.Parameters)` | cartesian `X`, `Y` | `z` argument accepted but ignored; propagation comes from `obj.Parameters` | yes, `Waist`, `PhiPhase`, and modal indices | needs-fix |
| `LaguerreBeam` | `opticalField(obj, X, Y, z)` -> `cart2pol(X, Y)` -> `computeComplexField(R, TH, obj.Parameters)` | cartesian `X`, `Y` converted to polar `R`, `TH` | `z` argument accepted but ignored; propagation comes from `obj.Parameters` | yes, `Waist`, `PhiPhase`, and modal indices | needs-fix |
| `ElegantHermiteBeam` | `opticalField(obj, X, Y, z)` -> `computeComplexField(X, Y, obj.Parameters)` | cartesian `X`, `Y` | `z` argument accepted but ignored; propagation comes from `obj.Parameters` | yes, `alpha`, `PhiPhase`, and modal indices | needs-fix |
| `ElegantLaguerreBeam` | `opticalField(obj, X, Y, z)` -> `cart2pol(X, Y)` -> `computeComplexField(R, TH, obj.Parameters)` | cartesian `X`, `Y` converted to polar `R`, `TH` | `z` argument accepted but ignored; propagation comes from `obj.Parameters` | yes, `alpha`, `PhiPhase`, and modal indices | needs-fix |

**Files:**
- Reference: `ParaxialBeams/ParaxialBeam.m`
- Reference: `ParaxialBeams/GaussianBeam.m`
- Reference: `ParaxialBeams/HermiteBeam.m`
- Reference: `ParaxialBeams/LaguerreBeam.m`
- Reference: `ParaxialBeams/ElegantHermiteBeam.m`
- Reference: `ParaxialBeams/ElegantLaguerreBeam.m`
- Modify: `README.md`
- Modify: `plan.md`

**Step 1: Write the failing audit checklist**

Agregar en `plan.md` una tabla `API Audit Matrix` con columnas `Class`, `Primary field method`, `Accepted coordinates`, `z semantics`, `Uses Parameters state`, `Status`.

**Step 2: Run the audit against the real classes**

Run: leer los seis archivos de beam y completar la tabla solo con evidencia del código.
Expected: toda fila queda marcada como `aligned`, `document-only`, o `needs-fix`.

**Step 3: Document the canonical contract**

Agregar en `README.md` una sección corta con este formato:

```md
## Beam API Contract

- canonical field entrypoint: `opticalField(...)`
- `Parameters` define constantes/modelo, no reemplaza argumentos dinámicos ambiguamente
- cada beam debe declarar sus coordenadas aceptadas
- cualquier desvío temporal queda documentado hasta el post-merge cleanup
```

**Step 4: Re-run the audit summary**

Run: verificar que `README.md` y `plan.md` usen el mismo contrato textual.
Expected: no hay contradicciones entre plan y documentación pública.

**Step 5: Commit**

```bash
git add plan.md README.md
git commit -m "docs: document public beam api contract"
```

### Task 3: Classify Canonical Examples

| File | Tier | Reason | Action |
| --- | --- | --- | --- |
| `examples/MainGauss_refactored.m` | canonical | usa `PhysicalConstants`, `GridUtils`, `FFTUtils` y flujo moderno explícito para Gaussian beam | mantener como punto de entrada principal para propagación escalar |
| `examples/MainMultiMode.m` | canonical | demuestra Hermite y Laguerre con clases modernas y setup corto, legible y actual | mantener como entrada canónica para modos múltiples |
| `ExampleRayTracing.m` | canonical | cubre la superficie moderna de ray tracing integrada con `GaussianBeam`, `RayBundle` y `RayTracer` | mantener como entrada canónica para ray tracing |
| `examples/MainGauss.m` | legacy-usable | nombre visible pero anterior al script refactorizado y menos alineado con utilidades nuevas | no presentarlo como default; dejarlo como referencia histórica utilizable |
| `examples/MainHermite.m` | historical | convive con múltiples variantes de tesis/paper y no está señalado como camino moderno principal | no enlazar en onboarding público |
| `examples/MainLaguerre.m` | historical | forma parte de una familia de scripts heredados con variantes específicas y naming inconsistente | no enlazar en onboarding público |

**Files:**
- Reference: `examples/MainGauss_refactored.m`
- Reference: `examples/MainMultiMode.m`
- Reference: `ExampleRayTracing.m`
- Modify: `README.md`
- Modify: `plan.md`

**Step 1: Create the example classification table**

Agregar en `plan.md` una tabla `Example Classification` con columnas `File`, `Tier`, `Reason`, `Action`.

**Step 2: Audit the candidate examples**

Run: leer los scripts candidatos y clasificarlos como `canonical`, `legacy-usable`, o `historical`.
Expected: hay entre 3 y 5 ejemplos `canonical`.

**Step 3: Point the public docs to canonical entrypoints only**

Actualizar `README.md` para que la sección de uso apunte solo a los ejemplos `canonical` elegidos.

**Step 4: Verify no legacy example is presented as the default**

Run: revisar `README.md`.
Expected: el primer camino de entrada del usuario nuevo pasa por ejemplos auditados, no por scripts históricos ambiguos.

**Step 5: Commit**

```bash
git add plan.md README.md
git commit -m "docs: classify canonical examples"
```

### Task 4: Validate Critical Test Coverage

## Critical Coverage Gates

- [ ] `z = 0` no produce NaN/Inf inválidos
- [ ] modo cero Hermite/Laguerre mantiene equivalencia razonable con Gaussian
- [ ] ray tracing cilíndrico sigue estable
- [ ] `tests/test_all.m` corre vía `portable_runner()`
- [ ] el runner funciona en Octave y MATLAB

Gap -> portable runner coverage -> pre-merge: `portable_runner()` hoy solo ejecuta `test_PhysicalConstants.m` y `test_RayTracing.m`; los checks de beams (`z = 0`, modos cero, elegant variants) existen como tests individuales pero no entran todavía por `tests/test_all.m`.

Gap -> MATLAB verification environment -> pre-merge: la máquina actual solo tiene MATLAB Runtime (`D:\Program Files\MATLAB\MATLAB Runtime\R2025b`), no `matlab.exe`; la compatibilidad MATLAB queda sin evidencia ejecutable local hasta correr `matlab -batch "run('tests/test_all.m')"` en un entorno con MATLAB completo.

**Files:**
- Reference: `tests/test_all.m`
- Reference: `tests/portable_runner.m`
- Reference: `tests/test_GaussianBeam.m`
- Reference: `tests/test_HermiteBeam.m`
- Reference: `tests/test_LaguerreBeam.m`
- Reference: `tests/test_ElegantHermiteBeam.m`
- Reference: `tests/test_ElegantLaguerreBeam.m`
- Reference: `tests/test_RayTracing.m`
- Modify: `plan.md`
- Modify: `tests/README.md`

**Step 1: Write the failing coverage checklist**

Agregar en `plan.md` una checklist `Critical Coverage Gates` con estos checks exactos:

```md
- [ ] `z = 0` no produce NaN/Inf inválidos
- [ ] modo cero Hermite/Laguerre mantiene equivalencia razonable con Gaussian
- [ ] ray tracing cilíndrico sigue estable
- [ ] `tests/test_all.m` corre vía `portable_runner()`
- [ ] el runner funciona en Octave y MATLAB
```

**Step 2: Run the portable suite in Octave**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: exit code `0` y mensaje `=== ÉXITO: Todos los tests pasaron ===`.

**Step 3: Run the suite in MATLAB**

Run: `matlab -batch "run('tests/test_all.m')"`
Expected: ejecución sin excepción y suite verde.

**Step 4: Document any remaining coverage gaps**

Si algún check de la checklist no está cubierto, registrar el gap en `plan.md` con formato `Gap -> owner -> pre-merge/post-merge`.

**Step 5: Align the testing docs**

Actualizar `tests/README.md` para que el comando principal coincida con el runner real y los coverage gates auditados.

**Step 6: Commit**

```bash
git add plan.md tests/README.md
git commit -m "test: lock critical pre-merge coverage gates"
```

### Task 5: Align Repository Narrative

## Stale Claims To Correct

- `README.md` describe clases como carpetas `@.../`, pero `ParaxialBeams/` hoy expone archivos `.m` planos
- `README.md` documenta `BeamSimulation`, pero ese archivo no existe en `ParaxialBeams/`
- `README.md` muestra scripts raíz `MainGauss.m`, `MainHermite.m`, `MainLaguerre.m`, pero los ejemplos auditados viven en `examples/` y el entrypoint recomendado para ray tracing está en `ExampleRayTracing.m`
- `tests/README.md` ya no debe sugerir que MATLAB Runtime alcanza para correr la suite

**Files:**
- Modify: `README.md`
- Modify: `tests/README.md`
- Modify: `plan.md`

**Step 1: Identify stale claims**

Run: revisar `README.md` y `tests/README.md` buscando estructura inexistente, ejemplos obsoletos o comandos de test inconsistentes.
Expected: lista explícita de claims a corregir.

**Step 2: Apply minimal narrative fixes**

Corregir únicamente estas superficies:

```md
- estructura real del repo
- clases realmente presentes en `ParaxialBeams/`
- ejemplos recomendados hoy
- compatibilidad MATLAB/Octave
- forma correcta de correr tests
```

**Step 3: Verify docs vs repository structure**

Run: comparar `README.md`, `tests/README.md`, `ParaxialBeams/` y `tests/`.
Expected: ningún archivo documenta rutas o clases inexistentes.

**Step 4: Commit**

```bash
git add README.md tests/README.md plan.md
git commit -m "docs: align repository narrative with current structure"
```

### Task 6: Prepare Git Integration Hygiene

## Merge Delta Summary

- architecture
- tests
- ci
- portability
- examples
- docs

## Merge Strategy Decision

- decision: `merge commit`
- rationale: `integration/pre-master` ya agrupa una secuencia visible de checkpoints de hardening (scope, API, examples, tests, docs) y preservarlos da mejor trazabilidad que un `squash` único antes del merge a `main/master`

## Integration Commit Draft

Merge `integration/pre-master`: stabilize modern beam API, portability, tests, and canonical examples

- consolidates modern `ParaxialBeams/` classes and portability fixes
- preserves current test surface and runner compatibility
- documents canonical examples and deferred post-merge redesign work

**Files:**
- Modify: `plan.md`

**Step 1: Capture the merge delta categories**

Agregar en `plan.md` una sección `Merge Delta Summary` con estas categorías:

```md
- architecture
- tests
- ci
- portability
- examples
- docs
```

**Step 2: Evaluate merge strategy**

Run: comparar claridad de `merge commit` vs `squash` para `integration/pre-master`.
Expected: decisión explícita y rationale corto en `plan.md`.

**Step 3: Draft the integration commit message**

Agregar en `plan.md` un borrador con esta estructura:

```md
Merge `integration/pre-master`: stabilize modern beam API, portability, tests, and canonical examples

- consolidates modern `ParaxialBeams/` classes and portability fixes
- preserves current test surface and runner compatibility
- documents canonical examples and deferred post-merge redesign work
```

**Step 4: Commit**

```bash
git add plan.md
git commit -m "docs: prepare merge integration checklist"
```

### Task 7: Final Readiness Gate

**Files:**
- Modify: `plan.md`
- Reference: `README.md`
- Reference: `tests/README.md`

**Step 1: Create the final acceptance checklist**

Agregar al final de `plan.md` esta checklist exacta:

```md
- [ ] merge scope frozen
- [ ] public API audited and documented
- [ ] canonical examples chosen
- [ ] critical tests green in Octave and MATLAB
- [ ] repo docs aligned with reality
- [ ] merge strategy decided
- [ ] deferred redesign work documented separately
```

**Step 2: Perform the final doc review**

Run: leer `plan.md`, `README.md` y `tests/README.md` de punta a punta.
Expected: los tres documentos cuentan la misma historia técnica.

**Step 3: Mark unresolved items explicitly**

Si algo no está listo para merge, dejarlo como unchecked y agregar nota `Blocker:` o `Deferred:` en `plan.md`.

**Step 4: Commit**

```bash
git add plan.md README.md tests/README.md
git commit -m "docs: add final pre-merge readiness gate"
```

### Task 8: Execute the Merge

**Files:**
- Reference: `plan.md`

**Step 1: Verify all readiness checks are complete**

Run: revisar la checklist final de `plan.md`.
Expected: todos los items pre-merge están completos o los blockers están explícitos.

**Step 2: Merge the integration branch**

Run: `git checkout master && git merge --no-ff integration/pre-master`
Expected: merge limpio o conflictos acotados y entendibles.

**Step 3: Verify the merged branch still passes the portable suite**

Run: `octave --no-gui --eval "run('tests/test_all.m')"`
Expected: suite verde también después del merge.

**Step 4: Commit or finalize merge message**

```bash
git commit
```

Usar el mensaje preparado en la Task 6 si Git abre editor por merge commit.

### Task 9: Open the Post-Merge Track

**Files:**
- Modify: `plan.md`

**Step 1: Create the deferred work list**

Agregar en `plan.md` una sección `Post-Merge Track` con estos items:

```md
- OO cleanup: separar modelo vs cálculo
- package migration a `+paraxial/...`
- legacy policy para ejemplos históricos
- application/services layer para simulaciones completas
```

**Step 2: Verify no deferred task leaked into pre-merge acceptance**

Run: releer `plan.md`.
Expected: el documento separa claramente merge readiness de evolución futura.

**Step 3: Commit**

```bash
git add plan.md
git commit -m "docs: define post-merge redesign track"
```

Plan complete and saved to `plan.md`. Two execution options:

**1. Subagent-Driven (this session)** - I dispatch fresh subagent per task, review between tasks, fast iteration

**2. Parallel Session (separate)** - Open new session with executing-plans, batch execution with checkpoints

**Which approach?**
