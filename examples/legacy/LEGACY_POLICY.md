# Legacy Examples Policy

## Overview

This document defines the policy for managing legacy examples in `examples/`.

## Classification

### Canonical Examples (`examples/canonical/`)

**Purpose**: Recommended entrypoints for new users.

| File | Description |
|------|-------------|
| `MainGauss_refactored.m` | Gaussian beam propagation with modern API |
| `MainMultiMode.m` | Multi-mode Hermite/Laguerre demonstration |
| `ExampleRayTracing.m` | Ray tracing visualization |
| `ExampleHankelPropagation.m` | Hankel beam propagation with unified API |

**Characteristics**:
- Use `BeamFactory.create()` for beam instantiation
- Use `setpaths()` for path initialization
- Runnable in Octave 11.1.0+ and MATLAB R2020b+
- Well-documented with inline comments
- Output verified against theory

### Legacy Folders

#### `examples/legacy/research/`

**Purpose**: Scripts used in thesis or research papers with specific parameters.

**Characteristics**:
- Hardcoded research-specific parameters
- May produce thesis/paper-specific visualizations
- Not intended for general use
- Preserved for reproducibility of published results

**Files**: Currently empty (no thesis-specific scripts identified)

#### `examples/legacy/generators/`

**Purpose**: Scripts that generate specific figures or data for papers.

**Characteristics**:
- Depend on pre-computed data or external inputs
- Produce publication-quality figures
- Not standalone runnable without specific data

| File | Description |
|------|-------------|
| `GenerateAbsHermiteHankels.m` | Generate abs Hermite-Hankel plots |
| `GenerateHankelsHermitePlots.m` | Generate Hankel-Hermite figures |
| `GenerateLastFigurePaper.m` | Paper figure generator |
| `GenerateLastFigurePaper2Obstruction.m` | Paper figure with obstruction |
| `GenerateLateralView.m` | Lateral view generator |
| `HankelHermiteSlices.m` | Hankel-Hermite slice visualization |
| `HankelLaguerrePropagation.m` | Hankel-Laguerre propagation figures |

#### `examples/legacy/archive/`

**Purpose**: Old API examples retained for reference.

**Characteristics**:
- Use old constructor patterns (`GaussianBeam(R, params)`)
- Use legacy `addpath` patterns
- Not recommended for new code
- May contain useful code patterns for reference

| File | Description |
|------|-------------|
| `MainGauss.m` | Old Gaussian example |
| `MainHermite.m` | Old Hermite example (825 lines) |
| `MainLaguerre.m` | Old Laguerre example |
| `MainLaguerre2.m` | Alternative Laguerre script |
| `MainLaguerreObstructionPropagation.m` | Laguerre with obstruction |
| `MainAnalyticPropagationGauss.m` | Analytic Gaussian propagation |
| `MainAnalyticPropagationHermite.m` | Analytic Hermite propagation |
| `MainAnalyticPropagationLaguerre.m` | Analytic Laguerre propagation |
| `MainAnalyticPropagationeHermite.m` | Analytic elegant Hermite |
| `MainHermite1DPlots.m` | 1D Hermite plots |
| `MainHermitePaper.m` | Paper-specific Hermite |
| `MainHermitePropagateHankel.m` | Hermite-Hankel propagation |
| `MainHermiteTEST.m` | Hermite test script |
| `PropagationTest.m` | Propagation test |
| `slices.m` | Slice visualization (heavy) |
| `SlicesGauss.m` | Gaussian slices |
| `SlicesHermite.m` | Hermite slices |
| `SlicesRays.m` | Ray slices |
| `AnalysisExperimentalLaguerre.m` | Experimental Laguerre analysis |
| `MaineLaguerreForTest.m` | Laguerre test variant |
| `MaineLaguerreWithOutObstruction.m` | Laguerre without obstruction |
| `MainLaguerreForTest.m` | Laguerre test variant |

## Migration Path

### Upgrading Legacy Scripts

To upgrade a legacy script to canonical:

```matlab
% BEFORE (old API)
addpath('ParaxialBeams')
addpath('ParaxialBeams', 'Addons')
GB = GaussianBeam(R, GaussianBeamParameters);

% AFTER (modern API)
setpaths()
GB = BeamFactory.create('gaussian', w0, lambda);
field = GB.opticalField(X, Y, z);
```

### New Script Submission

When adding new examples:

1. **Canonical**: Must use `BeamFactory.create()` and `setpaths()`
2. **Research**: Place in appropriate legacy folder with documentation
3. **Generators**: Place in `legacy/generators/` with data dependency notes

## Deprecation Policy

Scripts in `legacy/archive/` are:
- NOT actively maintained
- Kept for historical reference
- May be removed in future major versions

Scripts in `legacy/research/` and `legacy/generators/` are:
- Preserved for reproducibility
- NOT actively maintained
- May be archived to external storage in future

## References

- Migration plan: `docs/migration/LEGACY_MIGRATION_PLAN.md`
- Architecture: `docs/ARCHITECTURE.md`
- Package migration: `+paraxial/README.md`