# +paraxial Package (Future Development)

This directory contains the future MATLAB package structure for `paraxial`. 

## Status

**⚠️ This is a placeholder for future development. The current code uses individual `.m` files in `src/` and is NOT yet migrated to this package structure.**

## Planned Structure

```text
+paraxial/
├── +beams/                  % Beam classes
│   ├── ParaxialBeam.m      % Abstract base class
│   ├── GaussianBeam.m
│   ├── HermiteBeam.m
│   ├── LaguerreBeam.m
│   ├── ElegantHermiteBeam.m
│   ├── ElegantLaguerreBeam.m
│   ├── HankelHermite.m
│   └── HankelLaguerre.m
├── +parameters/             % Beam parameters
│   ├── GaussianParameters.m
│   ├── HermiteParameters.m
│   ├── LaguerreParameters.m
│   ├── ElegantHermiteParameters.m
│   └── ElegantLaguerreParameters.m
├── +propagation/
│   ├── +field/             % Field-based propagation
│   │   ├── IPropagator.m
│   │   ├── FFTPropagator.m
│   │   └── AnalyticPropagator.m
│   └── +rays/              % Ray-based propagation
│       ├── RayTracePropagator.m
│       ├── RayBundle.m
│       ├── RayTracer.m
│       ├── OpticalRay.m
│       └── CylindricalRay.m
└── +visualization/
    └── VisualizationUtils.m
```

## Migration Notes

### Why Package Migration?

Moving to `+paraxial/` package structure provides:
- **Namespace isolation**: `paraxial.GaussianBeam` vs `GaussianBeam`
- **Better collision avoidance**: No conflicts with other toolboxes
- **Cleaner imports**: `import paraxial.*` or `import paraxial.beams.*`

### Migration Path

1. Current state: Individual `.m` files in `src/`
2. Next step: Keep `src/` for compatibility, develop `+paraxial/` in parallel
3. Final step: Deprecate `src/` once `+paraxial/` is stable

### Usage Example (Future)

```matlab
% Import the entire package
import paraxial.*

% Create a beam using package namespace
beam = paraxial.beams.GaussianBeam(w0, lambda);

% Or import specific submodule
import paraxial.beams.GaussianBeam
beam = GaussianBeam(w0, lambda);
```

## Current vs Future

| Aspect | Current (`src/`) | Future (`+paraxial/`) |
|--------|------------------|------------------------|
| Namespace | Global | `paraxial.*` |
| Import | Manual addpath | `import paraxial.*` |
| Collision risk | High | Low |
| MATLAB best practice | Legacy | Modern |

## References

- [MATLAB Packages Documentation](https://mathworks.com/help/matlab/matlab_prog/create-and-use-packages.html)
- [MATLAB Packages Namespaces](https://mathworks.com/help/matlab/matlab_oop/package-overview.html)
