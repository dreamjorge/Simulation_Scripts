# Design: parameters-cleanup-sdd

## Technical Approach

Use a delegation-first refactor:

- Keep current parameter classes as stable public façades.
- Move numerical formulas to a dedicated stateless computation layer.
- Preserve all existing entrypoints while redirecting internal logic.

This minimizes risk and allows incremental verification.

## Architecture Decisions

### Decision: Delegation vs full DTO rewrite

**Choice**: Delegation (façade classes calling stateless utilities)

**Alternatives considered**:
- Full DTO-only parameter objects + API-breaking migration

**Rationale**:
- Required constraints include API stability and compatibility.
- Delegation yields separation of concerns without public breakage.

### Decision: New computation folder in `src/computation/`

**Choice**: Introduce dedicated folder for pure computation classes.

**Alternatives considered**:
- Keep formulas in `PhysicalConstants`
- Keep formulas spread across parameter classes

**Rationale**:
- Explicit architectural boundary improves discoverability and extension.
- Avoids further overloading `PhysicalConstants` semantics.

### Decision: Legacy helper extraction with shim

**Choice**: Move Hermite legacy helper implementation to `HermiteComputation` and keep compatibility shim in `HermiteParameters`.

**Alternatives considered**:
- Remove old entrypoint immediately

**Rationale**:
- Keeps backward compatibility for historical scripts while enforcing new ownership.

## Data Flow

### Current (before)

`Beams/Propagators -> Parameter classes (data + formulas + legacy helpers)`

### Target (after)

`Beams/Propagators -> Parameter classes (API façade/state) -> Computation layer (pure formulas)`

## Component Design

### `src/computation/BeamComputation.m`

Static methods (initial set):
- `rayleighDistance(w0, lambda)`
- `waveNumber(lambda)`
- `waist(w0, z, lambda, zr)`
- `gouyPhase(z, zr)`
- `radiusOfCurvature(z, zr)`
- `complexBeamParameter(z, zr, k)`

### `src/computation/HermiteComputation.m`

Static methods:
- `hermiteSolutions(nu, x)` (migrated legacy implementation)

### Parameter class updates

- `GaussianParameters`: delegate constructor-calculated and dynamic evaluations.
- `HermiteParameters`: delegate HG-related formulas and keep shim for `getHermiteSolutions`.
- `LaguerreParameters`: delegate LG formulas.
- `ElegantHermiteParameters` / `ElegantLaguerreParameters`: delegate alpha and modal phase formulas.

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/computation/BeamComputation.m` | Create | Centralized stateless formulas |
| `src/computation/HermiteComputation.m` | Create | Legacy Hermite helper owner |
| `src/parameters/GaussianParameters.m` | Modify | Delegation wiring |
| `src/parameters/HermiteParameters.m` | Modify | Delegation + shim |
| `src/parameters/LaguerreParameters.m` | Modify | Delegation wiring |
| `src/parameters/ElegantHermiteParameters.m` | Modify | Delegation wiring |
| `src/parameters/ElegantLaguerreParameters.m` | Modify | Delegation wiring |
| `setpaths.m` | Modify | Add computation path |
| `tests/portable_runner.m` | Modify | Add computation path |
| `tests/modern/test_BeamComputation.m` | Create | Core formula tests |
| `tests/modern/test_*Parameters.m` | Modify | Equivalence/compat tests |
| `tests/legacy_compat/*` | Modify | Shim compatibility checks |
| `docs/ARCHITECTURE.md` | Modify | New dependency direction |

## Testing Strategy

1. Focused formula tests for `BeamComputation`.
2. Parameter regression/equivalence tests (snapshot + dynamic API).
3. Full portable suite (`tests/test_all.m`).
4. Legacy compatibility suite (`tests/legacy_compat/run_legacy_compat.m`).

## Open Questions

- Whether to route `PhysicalConstants` internals to `BeamComputation` now or leave as separate static utility.
- Whether to add deprecation warnings for shim methods in this change or defer to next cleanup cycle.
