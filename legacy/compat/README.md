# Legacy Compatibility Layer

This directory contains deprecated aliases and adapters for backward compatibility with legacy scripts.

## Contents

### Deprecated Aliases

| Legacy Name | Modern Replacement | Notes |
|-------------|-------------------|-------|
| `HankeleHermite` | `HankelHermite` | Emits deprecation warning |
| `HankeleLaguerre` | `HankelLaguerre` | Emits deprecation warning |

## Usage

These aliases are retained for backward compatibility with historical scripts. **New code should use the modern classes directly.**

### MATLAB/Octave

```matlab
% Add modern library paths
addpath('src/beams', 'src/parameters', 'src/propagation/field', 'src/propagation/rays', 'src/visualization');

% Add legacy compatibility layer (optional - emits warnings)
addpath('legacy/compat');

% Modern approach (recommended)
h = HankelHermite(x, y, params, type);

% Legacy approach (deprecated, emits warning)
h = HankeleHermite(x, y, params, type);
```

## Migration Notes

See `docs/migration/LEGACY_MIGRATION_PLAN.md` for complete migration instructions.

## Deprecation Warnings

When using these aliases, you will see a warning:

```
warning: HankeleHermite is deprecated. Use HankelHermite instead.
warning: HankeleLaguerre is deprecated. Use HankelLaguerre instead.
```

This is expected behavior. The warnings help identify code that needs to be updated to the modern API.
