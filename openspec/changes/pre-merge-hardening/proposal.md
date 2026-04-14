# Proposal: Pre-Merge Hardening — `integration/pre-master` → `master`

## Intent

Consolidar y documentar el estado actual de `integration/pre-master` para un merge limpio a `master`. El branch tiene 71 commits de refactorización masiva (Strategy/Factory patterns, ~380 tests, CI/CD para MATLAB/Octave) pero carece de documentación que refleje la realidad del código. Sin este trabajo, el merge perpetúa la desalineación entre docs y código.

## Scope

### In Scope
- Actualizar `README.md` con estructura real del repo (27 archivos .m en ParaxialBeams)
- Crear `docs/ARCHITECTURE.md` con diagrama de clases y patrones aplicados
- Documentar API pública de beams (contrato `opticalField`, `getParameters`, `beamName`)
- Clasificar ejemplos como `canonical` vs `legacy` en `README.md`
- Completar checklist de merge readiness en `plan.md`
- Agregar CHANGELOG.md o conventional commits al merge
- Unificar branches redundantes (`exec/pre-merge-hardening` vs `integration/pre-master`)

### Out of Scope
- Package migration a `+paraxial/` namespaces (post-merge)
- Re-diseño OO de beams/propagators (post-merge)
- Reescritura total de ejemplos históricos (post-merge)
- Rescate de ramas legacy completas

## Capabilities

### New Capabilities
- `pre-merge-docs`: Documentación de arquitectura y API para facilitar merge a master

### Modified Capabilities
- `beam-api`: El contrato existente de ParaxialBeam se documenta formalmente pero no cambia

## Approach

1. **Auditoría de API**: Leer los 6 archivos de beam classes y verificar que cumplen el contrato `opticalField(X,Y,z)`, `getParameters(z)`, `beamName()`
2. **Documentación canónica**: Crear `docs/ARCHITECTURE.md` con estructura de clases, patrones, y flujo de datos
3. **README.md rewrite**: Reemplazar estructura vieja ( `@Carpeta/` ) con realidad actual ( `.m` files)
4. **Clasificación de ejemplos**: Marcar `examples/MainGauss_refactored.m`, `examples/MainMultiMode.m`, `ExampleRayTracing.m` como canonical
5. **Unificación de branches**: Comparar `exec/pre-merge-hardening` con `integration/pre-master` y consolidar o archivar

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `README.md` | Modified | Estructura actual, ejemplos canonicales, MATLAB/Octave compatibility |
| `plan.md` | Modified | Completar readiness checklist pre-merge |
| `docs/ARCHITECTURE.md` | New | Diagrama de clases, patrones Strategy/Factory, flujo de datos |
| `examples/` | Modified | Clasificar scripts como canonical vs legacy |
| `openspec/` | New | Estructura SDD para tracking |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Conflictos entre `exec/pre-merge-hardening` e `integration/pre-master` | Medium | Hacer diff de ambas vs master y decidir cuál absorber |
| `README.md` desactualizado confunde usuarios post-merge | High | Verificación cruzada con `ls ParaxialBeams/*.m` antes de commit |
| Tests quebrados en MATLAB (solo verificado en Octave) | Low | CI ya tiene workflow de MATLAB |

## Rollback Plan

```bash
# Si el merge tiene problemas:
git revert <merge-commit>
git checkout master
# Restaurar docs desde el commit anterior al merge
```

## Dependencies

- GitHub Actions CI passing para Octave y MATLAB (ya configurado)
- 71 commits de refactor en `integration/pre-master` verificados

## Success Criteria

- [ ] `README.md` refleja estructura real: `ls ParaxialBeams/*.m | wc -l` = 27 archivos documentados
- [ ] `docs/ARCHITECTURE.md` existe con diagrama de arquitectura
- [ ] `plan.md` tiene todos los items de pre-merge checklist marcados ✅
- [ ] Ejemplos canonicales identificados y documentados
- [ ] `exec/pre-merge-hardening` archivado o mergeado
- [ ] Merge commit sigue conventional commits: `merge(integration): stabilize beam API...`
