# Design: Pre-Merge Hardening Documentation

## Technical Approach

Update and create documentation artifacts to accurately reflect the current `integration/pre-master` codebase before merging to `master`. This is documentation-only work—no code changes.

## Architecture Decisions

### Decision: README.md Rewrite vs Partial Update

**Choice**: Complete rewrite of README.md structure section
**Alternatives considered**: Incremental patches to existing README
**Rationale**: Existing README references `@Folder/` class folders that no longer exist. A complete rewrite prevents confusion from mixing old and new structures.

### Decision: Separate ARCHITECTURE.md vs Inline README Section

**Choice**: Create separate `docs/ARCHITECTURE.md`
**Alternatives considered**: Add architecture section to README.md
**Rationale**: Architecture is complex enough (6 beam types + 3 propagators + factory + strategy pattern) that a separate document with diagrams is clearer. README focuses on "how to use", ARCHITECTURE focuses on "how it works".

### Decision: Example Classification Criteria

**Choice**: Classify based on (1) executes without error, (2) documented as canonical, (3) maintained in sync with API
**Alternatives considered**: Classification by age or file size
**Rationale**: User trust depends on examples actually working. Age/size are poor proxies for reliability.

## Data Flow

Current system architecture:

```
┌─────────────────────────────────────────────────────────────────┐
│                        User Code                                │
│  examples/MainGauss_refactored.m, MainMultiMode.m, etc.        │
└─────────────────────────┬───────────────────────────────────────┘
                          │ addpath ParaxialBeams
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                      BeamFactory                                │
│  BeamFactory.create('gaussian', w0, lambda)                    │
└─────────────────────────┬───────────────────────────────────────┘
                          │ creates
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                   ParaxialBeam (abstract)                      │
│  opticalField(X,Y,z) │ getParameters(z) │ beamName()           │
└───────┬───────┬───────┬───────┬───────┬───────┬─────────────────┘
        │       │       │       │       │       │
   ┌────┴┐ ┌───┴───┐ ┌─┴──┐ ┌──┴──┐ ┌──┴──┐ ┌───┴────┐
   │Gauss│ │Hermit│ │Lag │ │ElegH│ │ElegL│ │HankelL │
   └─────┘ └──────┘ └────┘ └─────┘ └─────┘ └────────┘
        │       │       │       │       │       │
        └───────┴───────┴───────┴───────┴───────┘
                          │ propagate(beam, z)
                          ▼
┌─────────────────────────────────────────────────────────────────┐
│                   IPropagator (Strategy)                        │
│  ┌──────────────┐ ┌───────────────┐ ┌────────────────┐        │
│  │FFTPropagator│ │AnalyticPropag│ │RayTracePropagat│        │
│  └──────────────┘ └───────────────┘ └────────────────┘        │
└─────────────────────────────────────────────────────────────────┘
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `README.md` | Modify | Rewrite structure section, add canonical examples, update MATLAB/Octave compatibility |
| `docs/ARCHITECTURE.md` | Create | Class hierarchy diagram, pattern descriptions, data flow |
| `examples/*.m` | Modify | Add header comment: `%% canonical` or `%% legacy` |
| `plan.md` | Modify | Complete pre-merge checklist items |

## Testing Strategy

This is documentation-only. Verification by:
1. `ls ParaxialBeams/*.m | wc -l` confirms 27 files match README
2. Running canonical examples: `octave examples/MainGauss_refactored.m`
3. Reading `help BeamFactory` confirms documented API works

## Open Questions

- [ ] Should `exec/pre-merge-hardening` branch be archived or merged?
- [ ] Which examples beyond `MainGauss_refactored.m`, `MainMultiMode.m`, `ExampleRayTracing.m` qualify as canonical?
