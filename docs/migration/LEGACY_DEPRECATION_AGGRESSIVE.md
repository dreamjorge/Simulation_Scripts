# Legacy Alias Deprecation Plan (Agresivo)

**Status:** PLANNED
**Timeline:** 1 mes post-merge
**Affected:** `legacy/compat/HankeleHermite.m`, `legacy/compat/HankeleLaguerre.m`

---

## Background

Los aliases legacy (`HankeleHermite`, `HankeleLaguerre`) emiten warnings de migración pero aún funcionan. El plan de modernización requiere remoção agresiva post-merge para evitar technical debt.

---

## Deprecation Timeline

### Week 1 (Post-Merge Day 1-7): Warning → Error

**Cambios:**
- Modificar `legacy/compat/HankeleHermite.m` para que tire **error** en lugar de warning
- Modificar `legacy/compat/HankeleLaguerre.m` para que tire **error** en lugar de warning

**Código cambio en ambos archivos:**
```matlab
% ANTES (warning):
persistent warnIssued
if isempty(warnIssued)
    warning('HankeleHermite is deprecated...');
    warnIssued = true;
end

% DESPUÉS (error):
error('HankeleHermite has been removed. Use HankelHermite instead. See docs/migration/LEGACY_MIGRATION_PLAN.md');
```

**Verificación:**
```bash
git diff legacy/compat/HankeleHermite.m legacy/compat/HankeleLaguerre.m
octave --no-gui --eval "addpath('legacy/compat'); HankeleHermite()"
# Debe producir error (no warning)
```

---

### Week 2 (Post-Merge Day 8-14): Remove Files

**Cambios:**
- Eliminar `legacy/compat/HankeleHermite.m`
- Eliminar `legacy/compat/HankeleLaguerre.m`

**Verificación:**
```bash
ls legacy/compat/
# HankeleHermite.m y HankeleLaguerre.m NO deben existir
git status
# Debe mostrar deletion
```

---

### Week 3 (Post-Merge Day 15-21): Documentation Update

**Cambios:**
- Actualizar `docs/migration/LEGACY_MIGRATION_PLAN.md` para marcar aliases como **REMOVED**
- Actualizar `README.md` si menciona los aliases
- Agregar nota en `CHANGELOG.md`

---

## Risk Analysis

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Usuario existente tiene código con HankeleHermite | Medium | High (breakage) | Comunicar 1 mes antes via CHANGELOG y GitHub release notes |
| Tests en `legacy_compat/` fallan | High | Medium (CI rojo) | Remover tests de legacy_compat en Week 2 también |
| Dependency de otras classes en legacy/compat | Low | Medium | Verificar antes de remover |

---

## Affected Files

### Aliases to Remove
- `legacy/compat/HankeleHermite.m`
- `legacy/compat/HankeleLaguerre.m`

### Tests to Remove
- `tests/legacy_compat/test_HankelCompatibility.m` (si solo testa aliases)
- `tests/legacy_compat/test_HankelAliasEdgeCases.m`
- `tests/legacy_compat/test_HankelAliasStaticDelegation.m`

### Documentation to Update
- `docs/migration/LEGACY_MIGRATION_PLAN.md`
- `README.md` (si menciona aliases)
- `CHANGELOG.md`

---

## Communication Template

Para GitHub release notes / commit message:

```markdown
## Breaking Change: Legacy Aliases Removed

Los aliases `HankeleHermite` y `HankeleLaguerre` han sido eliminados.

**Acción requerida:** Usar `HankelHermite` y `HankelLaguerre` directamente.

Ver docs/migration/LEGACY_MIGRATION_PLAN.md para guía de migración.
```

---

## Verification Commands

```bash
# Week 1: Verify error is thrown
cd D:\Repositories\Simulation_Scripts
octave --no-gui --eval "addpath('legacy/compat'); HankeleHermite()"
# Expected: error message

# Week 2: Verify file removed
ls legacy/compat/
# Expected: HankeleHermite.m not found

# Week 3: Verify CI still green (sin legacy alias tests)
octave --no-gui --eval "addpath('tests'); portable_runner();"
# Expected: all tests pass (without legacy alias tests)
```

---

## Rollback Plan (si es necesario)

Si hay breakage crítico inesperado:
```bash
git checkout HEAD~1 -- legacy/compat/HankeleHermite.m legacy/compat/HankeleLaguerre.m
```

---

**Fecha de inicio:** Post-merge a master
**Duración total:** 3 semanas
