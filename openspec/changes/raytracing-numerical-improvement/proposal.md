# Proposal: raytracing-numerical-improvement

## Intent

Reemplazar el gradiente de fase por diferencias forward con delta fijo en `RayTracer` y `HankelRayTracer` por un método basado en el campo complejo que evita `unwrap` y mejora la precisión numérica. Simultáneamente corregir la inconsistencia de estado en `RayBundle` donde `r` y `theta` no se actualizan después de `addStep`, y arreglar el criterio frágil de axis-crossing en `HankelRayTracer`.

## Scope

### In Scope
- Reemplazar `calculateSlopes()` en `RayTracer.m` con gradiente complejo: `∇φ = Im(u̅∇u) / (|u|² + ε)`
- Agregar `resolveDelta()` privado escalado por waist, λ y escala local
- Convertir `r` y `theta` en `Dependent` en `RayBundle.m`
- Reemplazar criterio `(x0*y1 - x1*y0) < 0` en `HankelRayTracer` por distancia mínima del segmento al origen
- Agregar protección contra `|field|` pequeño en zonas de baja amplitud
- Agregar tests de exactitud física en `test_RayTracing.m`

### Out of Scope
- API analítica `phaseGradient()` en beams (futura capacidad)
- Política adaptativa de `dz` en `RayTracePropagator`
- Rework de `HankelLaguerre` o `HankelHermite` con API propia

## Capabilities

### New Capabilities
- `complex-phase-gradient`: Gradiente de fase calculado desde campo complejo sin `unwrap`, más estable cerca de singularidades

### Modified Capabilities
- `ray-slope-calculation`: Mejora de método numérico — precisión de gradiente de orden superior y escalado adaptativo del delta

## Approach

Reemplazar la diferencia forward hardcodeada por gradiente complejo regularizado. La identidad `∇φ = Im(u̅∇u)/|u|²` se implementa con diferencia central para el gradiente espacial y un ε regularizador para evitar división por cero en regiones de baja amplitud. El δ se calcula dinámicamente como `max(λ, |x| * 1e-4)` para adaptarse a la escala espacial local.

Para `RayBundle`, convertir `r` y `theta` en Dependent que 计算 desde `x` e `y` en cada acceso, eliminando el riesgo de estado inconsistente.

Para el axis-crossing, calcular la distancia mínima del segmento de rayo `(x0,y0)→(x1,y1)` al origen y verificar que la proyección del origen sobre el segmento caiga dentro del tramo.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `src/propagation/rays/RayTracer.m` | Modified | Nuevo método `calculatePhaseGradientComplex()`, `resolveDelta()` privado, protección ε |
| `src/propagation/rays/RayBundle.m` | Modified | `r` y `theta` pasan a Dependent, `addStep()` ya no los calcula |
| `src/propagation/rays/HankelRayTracer.m` | Modified | Criterio geométrico real de axis-crossing, mismo gradiente complejo |
| `tests/modern/test_RayTracing.m` | Modified | Tests de convergencia numérica vs analítico para Gaussian, verificación Euler vs RK4 |
| `tests/edge_cases/test_RayTracing_extreme.m` | Modified | Corregir tests que setearon `bundle.sx/sy` sin efecto |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Breaking existing beam consumers | Low | La firma de `calculateSlopes` no cambia, solo el algoritmo interno |
| Edge cases con | Med | ε regularizador y fallback a diferencia central si `|u| < ε` |

## Rollback Plan

Revertir los tres archivos `.m` a su estado anterior desde git. Los tests no son código de producción así que no requieren rollback separados.

## Dependencies

Ninguna dependencia externa. Todo es contenido dentro del repo.

## Success Criteria

- [ ] Tests de `test_RayTracing.m` siguen pasando con precisión numérica mejorada
- [ ] Tests de `test_HankelRayTracing.m` verifican flip correcto en axis-crossing real
- [ ] `test_RayTracing_extreme.m` no falla por cambios en cómo se setean slopes (reconciliar test vs implementación)
- [ ] Gradiente numérico converge al analítico para GaussianBeam cuando se refina δ
- [ ] `bundle.r` y `bundle.theta` siempre reflejan la última posición (no hay estado inconsistente)