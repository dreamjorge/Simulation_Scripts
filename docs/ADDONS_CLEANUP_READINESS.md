# Addons Cleanup Readiness

This document summarizes disposal/retention decisions for `ParaxialBeams/Addons/`
as a result of the deep-inventory phase. It is **not** permission to delete files
immediately — every removal requires a dedicated SDD/OpenSpec change with guardrail
coverage and passing tests.

---

## Keep During Compatibility Window

The following addons are referenced by legacy archive scripts, research examples,
or internal calls within the addon family itself. They must remain until the legacy
compatibility window is formally closed via an SDD.

| Addon | Rationale |
|-------|-----------|
| `getPropagateCylindricalRays.m` | Called by legacy archive scripts; internally calls `getCylindricalGradient`, `copyArrayRay2Ray`, `copyRay2ArrayRay`. Hankel class static API (`HankelLaguerre.getPropagateCylindricalRays`) provides modern equivalent. |
| `getCylindricalGradient.m` | Internal dependency of `getPropagateCylindricalRays.m`. Referenced in `AnalysisUtils.m` comments as matching legacy logic. |
| `assignCoordinates2CartesianRay.m` | Referenced by 20+ research/archive scripts and `OpticalRay` class comments. |
| `assignCoordinates2CylindricalRay.m` | Same as above — active use in legacy Laguerre archive scripts. |
| `copyRay2ArrayRay.m` | Referenced in `OpticalRay` comments as expected helper. Called internally by `getPropagateCylindricalRays.m`. |
| `copyArrayRay2Ray.m` | Same as above. |
| `copyRay.m` | Referenced by historical archive scripts; part of legacy copy helper family. |
| `copy2Ray.m` | Same family. Usage not proven outside legacy scripts, but kept pending investigation. |
| `copyElementRay.m` | Same family. Kept pending broader history search. |
| `copyElementsOnRay.m` | Same family. |
| `getPropagateRay.asv` | MATLAB autosave artifact. Must be removed only in a dedicated cleanup change. |

---

## Keep for Legacy Plot Reproducibility

These addons are used exclusively by `examples/legacy/` research and archive scripts
to produce published figures. Removing them would break historical reproducibility.

| Addon | Rationale |
|-------|-----------|
| `AdvancedColormap.m` | Referenced across research scripts (`MaineHermiteForThesisParameters.m`, `MainHermiteThesis.m`, generator scripts). Used for figure colormap generation. |
| `tight_subplot.m` | Same research scripts use it for multi-panel figure layout. |
| `paraxialPropagator.m` | Research scripts use for propagation loops. |
| `propagateOpticalField.m` | Same research scripts for field propagation. |
| `unwrap_phase.m` | Used by phase visualization scripts. |
| `Plots_Functions/` dir | Contains `plotOpticalField`, `plotGaussianParameters`, and ray plotting helpers for legacy examples. |

---

## License and Origin Review Required

These addons appear to be third-party vendored distributions. Before any relocation,
redistribution, or replacement decision, origin and license must be confirmed.

| Addon | Origin Clues | License Clues | Recommended Action |
|-------|--------------|---------------|-------------------|
| `export_fig-master/` | Unmodified vendored MATLAB export utility. No local attribution header. Header says "Downloaded from MathWorks File Exchange" (comment in `export_fig.m`). | Internal `license.txt` covering this tree. BSD-3 from John D'Errico. | Confirm FEX submission number and version. Decide to replace with `print`/`saveas` native alternatives or keep for legacy scripts. |
| `panel-2.14/` | Unmodified vendored panel layout utility. No attribution header. | `license.txt` applies here too (BSD-3 John D'Errico). | Consider native `uipanel` replacement for new code; keep for legacy compatibility. |
| `vline.m` | Common plotting utility. No attribution header. Likely File Exchange origin. | Unknown. No license file present in this directory. | Search FEX history for origin. Confirm no GPL/copyleft contamination before using in new code. |
| `tight_subplot.m` | Common File Exchange utility name. No attribution header. | Unknown. | Same as `vline.m`. |
| `unwrap_phase.m` | Common signal processing utility. No attribution header. | Unknown. | Verify origin before using in new code. |
| `license.txt` | BSD-3 John D'Errico 2012. Applies to `export_fig-master/` and `panel-2.14/`. | Explicit BSD-3 — permissive. | Keep with vendored assets. |

**Action item:** Create a dedicated cleanup SDD to address `vline`, `tight_subplot`, and `unwrap_phase` origin before using them in new code.

---

## Needs Usage Investigation

These addons were flagged in prior passes and require broader history analysis
before any reclassification.

| Addon | Rationale |
|-------|-----------|
| `copy2Ray.m` | Part of legacy copy helper family. Not directly referenced in `+paraxial/`, `src/`, or `tests/` — only in legacy archive scripts. |
| `copyElementRay.m` | Same family. External usage not proven in this pass. |
| `copyElementsOnRay.m` | Same family. |

---

## Potential Future Removal Candidates

The following are `removable-candidate` only after a formal SDD proves no active
references remain AND a guardrail test confirms no regressions.

- `paraxialPropagator.m` — after legacy examples are refactored to use `+paraxial/+propagation/+field/` equivalents.
- `propagateOpticalField.m` — after field propagation examples are migrated.
- `AdvancedColormap.m`, `tight_subplot.m`, `vline.m`, `unwrap_phase.m` — after origin confirmed and examples migrated to native alternatives.

**Removal gates:**
1. SDD/OpenSpec change created with explicit rationale.
2. Guardrail test asserts no references in `+paraxial/`, `src/`, `tests/`, or `examples/` (excluding `examples/legacy/archive/` which may retain references during transition).
3. All tests pass.
4. Legacy compatibility window formally closed.

---

## Required Follow-up SDDs

The following cleanup decisions need their own SDD/OpenSpec changes:

| SDD | Scope |
|-----|-------|
| `addons-vendored-licensing` | Review `vline`, `tight_subplot`, `unwrap_phase` origin and license. Decide whether to keep, replace, or archive. |
| `addons-legacy-plot-migration` | Migrate research figure scripts from vendored helpers to native alternatives. |
| `addons-copy-helpers-investigation` | Confirm `copy2Ray`, `copyElementRay`, `copyElementsOnRay` are dead code or still referenced. |
| `addons-autosave-cleanup` | Remove `getPropagateRay.asv` in a dedicated cleanup change. |

---

## Structural Invariants (do not break)

- `ParaxialBeams/Addons/` directory must not be deleted while `examples/legacy/` scripts reference its contents.
- `license.txt` must remain in place while vendored distributions (`export_fig-master/`, `panel-2.14/`) are present.
- `getPropagateRay.asv` must not be committed as a new file — existing artifact may be cleaned up in `addons-autosave-cleanup` SDD.