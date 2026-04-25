# Design: +paraxial/ Package Migration (Phase 3)

## Technical Approach

Implement Strangler Fig migration using BeamFactory as the routing mechanism. BeamFactory receives `+paraxial/` priority via `which()` lookup, falls back to `src/` with warning, and emits deprecation guidance when `src/` classes are used directly.

## Architecture Decisions

### Decision: BeamFactory as migration routing point

**Choice**: BeamFactory.create() resolves classes at runtime using `which()` lookup
**Alternatives considered**: Compile-time class path constants; config flag to switch namespaces
**Rationale**: `which()` adapts to partial migrations — if only GaussianBeam is migrated, BeamFactory routes correctly. No config file to maintain.

### Decision: Deprecation via warning(), not @decorator

**Choice**: `warning()` call in src/ class constructors
**Alternatives considered**: `@deprecation` MATLAB attribute (R2019a+ only); custom deprecation class
**Rationale**: `warning()` works in both Octave and MATLAB. The message guides users to both BeamFactory and direct `+paraxial/` path.

### Decision: src/ classes stay functional during transition

**Choice**: src/ classes emit deprecation warning but still instantiate normally
**Alternatives considered**: src/ classes throw error forcing immediate migration; src/ classes return empty/placeholder
**Rationale**: Backward compatibility — existing user code using direct `src/beams/GaussianBeam()` still works while they migrate.

## Data Flow

```
User: BeamFactory.create('gaussian', w0, lambda)
         │
         ▼
BeamFactory.create()
  ├── which('GaussianBeam') → finds +paraxial/+beams/GaussianBeam.m
  │      └── Instantiate from +paraxial/ ✅
  └── if not found → which('src/beams/GaussianBeam.m')
         └── warning('src/beams is deprecated...')
         └── Instantiate from src/ ✅
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `ParaxialBeams/BeamFactory.m` | Modify | Add `which()`-based routing; fallback with warning |
| `src/beams/GaussianBeam.m` | Modify | Add deprecation warning in constructor |
| `src/beams/HermiteBeam.m` | Modify | Add deprecation warning in constructor |
| `src/beams/LaguerreBeam.m` | Modify | Add deprecation warning in constructor |
| `src/beams/ElegantHermiteBeam.m` | Modify | Add deprecation warning in constructor |
| `src/beams/ElegantLaguerreBeam.m` | Modify | Add deprecation warning in constructor |
| `src/beams/HankelLaguerre.m` | Modify | Add deprecation warning in constructor |
| `src/beams/HankelHermite.m` | Modify | Add deprecation warning in constructor |
| `README.md` | Modify | Document `+paraxial/` as canonical namespace |

## Implementation Details

### BeamFactory routing logic (new create() method):

```matlab
function beam = create(type, w0, lambda, varargin)
    % Parse n, m, l, p, htype from varargin...

    switch lower(type)
    case 'gaussian'
        className = 'GaussianBeam';
        canonical = 'paraxial.beams.GaussianBeam';
        legacy = 'src/beams/GaussianBeam.m';
    % ... same for other types
    end

    % Try +paraxial/ first
    if exist(canonical, 'class')
        beam = feval(canonical, w0, lambda, varargin{:});
    % Fallback to src/ with warning
    elseif exist(legacy, 'file')
        warning('BeamFactory:deprecatedSrc', ...
            ['src/beams/%s is deprecated. ' ...
             'Use BeamFactory.create(''%s'', ...) or %s directly.'], ...
            className, type, canonical);
        beam = feval(className, w0, lambda, varargin{:});
    else
        error('BeamFactory:classNotFound', ...
            'Neither +paraxial nor src/ version of %s found.', className);
    end
end
```

### Deprecation warning in src/ constructors:

```matlab
function obj = GaussianBeam(arg1, arg2, varargin)
    obj = obj@ParaxialBeam();
    if nargin == 0, return; end
    % Emit warning only when actually instantiated
    warning('BeamFactory:deprecated', ...
        'src/beams/GaussianBeam is deprecated. Use BeamFactory.create(''gaussian'', ...) or +paraxial/+beams/GaussianBeam directly.');
    % ... rest of constructor
end
```

### Migration sequence:

```
1. GaussianBeam → verify tests → commit
2. HermiteBeam  → verify tests → commit
3. LaguerreBeam → verify tests → commit
4. ElegantHermiteBeam → verify tests → commit
5. ElegantLaguerreBeam → verify tests → commit
6. HankelLaguerre → verify tests → commit
7. HankelHermite → verify tests → commit
8. Deprecate src/beams/ in README
```

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Unit | BeamFactory routing logic | Mock `exist()` to test both paths |
| Integration | Full `BeamFactory.create()` for each type | `test_BeamFactory.m` already exists |
| Migration | `src/beams/` classes emit warning | Verify warning() called in constructor tests |
| Standalone | `+paraxial/` classes without `src/` on path | Test with isolated path |

## Migration / Rollout

**No data migration required.** This is pure code reorganization.

Rollout: Per-class migration with test verification after each. If any class fails tests, halt migration and investigate before continuing.

## Open Questions

- None — design follows directly from proposal and spec.