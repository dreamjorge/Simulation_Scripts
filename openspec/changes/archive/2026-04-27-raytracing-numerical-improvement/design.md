# Design: raytracing-numerical-improvement

## Technical Approach

Reemplazar el gradiente de fase por diferencias forward con `delta` fijo y `unwrap(angle(field))` por un gradiente complejo regularizado `∇φ = Im(u̅∇u) / (|u|² + ε)` implementado con diferencia central y `δ` escalado dinámicamente. Los cambios se concentran en `RayTracer.m`, `RayBundle.m` y `HankelRayTracer.m` siguiendo el orden de prioridad del documento de análisis técnico.

## Architecture Decisions

### Decision: Complex gradient over forward difference with unwrap

**Choice**: Gradiente complejo `Im(u̅ ∇u) / |u|²` con ε regularización
**Alternatives considered**:
- Forward difference + `unwrap` (actual) — frágil cerca de singularidades
- Central difference con `unwrap` — mejor orden pero sigue dependiendo de unwrap
**Rationale**: Elimina dependencia de unwrap y es matemáticamente más cercano a la física del campo complejo. El ε maneja división por cero.

### Decision: Delta scaling formula

**Choice**: `δ = max(λ, |x|·1e-4, |y|·1e-4, w0·1e-4)`
**Alternatives considered**:
- Fixed `δ = 1e-7` (actual) — no adapta a escala
- `δ = λ/100` — ignora escala geométrica local
**Rationale**: Combina escala de longitud de onda con escala espacial local y waist. Suficientemente chico para curvatura fina, grande para evitar cancelación numérica.

### Decision: r/theta as Dependent vs recomputing in addStep

**Choice**: Dependent properties computed from x, y on access
**Alternatives considered**:
- Recalcular en `addStep()` — frágil, puede desincronizar
- Mantener como stored — requiere reescribir constructor
**Rationale**: Elimina clase entera de bugs de estado. MATLAB/Octave recalculan transparentemente.

### Decision: Axis crossing via minimum segment distance

**Choice**: Distancia mínima del segmento `(x0,y0)→(x1,y1)` al origen + proyección dentro del tramo
**Alternatives considered**:
- Signo de `(x0·y1 - x1·y0)` (actual) — solo mide orientación
- Distancia euclidiana a origen — falla en ángulos rasantes
**Rationale**: Geométricamente correcto. Detecta cruce real, no solo cambio de orientación.

## Data Flow

```
RayTracer.propagate()
│
├── get current (x0, y0, z0) from bundle
├── calculateSlopes(beam, x0, y0, z0)
│   │
│   ├── resolveDelta(x0, y0, w0, λ) → δ
│   ├── field = beam.opticalField(x0, y0, z0)
│   ├── field_dx = beam.opticalField(x0+δ, y0, z0)
│   ├── field_dy = beam.opticalField(x0, y0+δ, z0)
│   ├── ∂u/∂x ≈ (field_dx - field_dx') / (2δ)  [central diff]
│   ├── ∂u/∂y ≈ (field_dy - field_dy') / (2δ)
│   ├── sx = Im(conj(field) · ∂u/∂x) / (|field|² + ε)
│   └── sy = Im(conj(field) · ∂u/∂y) / (|field|² + ε)
├── RK4: k1,k2,k3,k4 stages
└── bundle.addStep(x1, y1, z1, sx, sy)
```

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/propagation/rays/RayTracer.m` | Modify | Agregar `calculatePhaseGradientComplex()`, `resolveDelta()`; modificar `calculateSlopes()` para usar complejo |
| `src/propagation/rays/RayBundle.m` | Modify | `r` y `theta` → Dependent con getter; `addStep()` ya no los calcula |
| `src/propagation/rays/HankelRayTracer.m` | Modify | Reemplazar criterio axis-crossing con distancia mínima; usar mismo gradiente complejo |
| `tests/modern/test_RayTracing.m` | Modify | Agregar tests convergencia vs analítico, Euler vs RK4 con tolerancias |
| `tests/edge_cases/test_RayTracing_extreme.m` | Modify | Corregir tests que setearon sx/sy sin efecto real |

## Interfaces / Contracts

### RayTracer.calculatePhaseGradientComplex (new private method)

```matlab
function [sx, sy] = calculatePhaseGradientComplex(beam, x, y, z)
    % Compute phase gradient: ∇φ = Im(u̅∇u) / (|u|² + ε)
    % Uses central difference for ∂u/∂x, ∂u/∂y
    % Returns sx, sy with same units as dx/dz, dy/dz
end
```

### RayTracer.resolveDelta (new private static method)

```matlab
function delta = resolveDelta(x, y, w0, lambda)
    % Adaptive delta: max(lambda, |x|*1e-4, |y|*1e-4, w0*1e-4)
    delta = max(lambda, abs(x) * 1e-4, abs(y) * 1e-4, w0 * 1e-4);
end
```

### RayBundle.r, theta (changed to Dependent)

```matlab
properties (Dependent)
    r
    theta
end

methods
    function val = get.r(obj)
        val = sqrt(obj.x(:,:,end).^2 + obj.y(:,:,end).^2);
    end
    function val = get.theta(obj)
        [~, val] = cart2pol(obj.x(:,:,end), obj.y(:,:,end));
    end
end
```

### HankelRayTracer axis crossing check (replaces line 60)

```matlab
% Distance from segment to origin
dx = x1 - x0; dy = y1 - y0;
t = -(x0.*dx + y0.*dy) ./ (dx.^2 + dy.^2 + eps);
t_clamped = max(0, min(1, t));
nearest_x = x0 + t_clamped .* dx;
nearest_y = y0 + t_clamped .* dy;
min_dist = sqrt(nearest_x.^2 + nearest_y.^2);
crossed = min_dist < max(abs(x0), abs(y0)) * 1e-3;
```

## Testing Strategy

| Layer | What to Test | Approach |
|-------|-------------|----------|
| Unit | `calculatePhaseGradientComplex` accuracy | Comparar vs analítico para Gaussian: `dφ/dx = k·x/R(z)`已知解 |
| Unit | `resolveDelta` scaling | Verificar δ ∝ λ cerca de origen, δ ∝ |x| lejos |
| Unit | Axis crossing geometric | Crear rayos que cruzan cerca del origen, verificar flip en momento correcto |
| Integration | `RayTracer.propagate` convergence | Refinar dz, verificar que Euler→RK4 converge |
| Integration | `RayBundle` state consistency | Acceder r/theta post-addStep, comparar con recalculado manual |

## Migration / Rollback

No migration required. Changes son contained a archivos existentes. Rollback vía git checkout de los tres archivos `.m`.

## Open Questions

- [ ] ¿Debería `resolveDelta` usar también `dy` (no solo `y`) para grids no simétricos?
- [ ] ¿El ε regularizador debe ser absoluto (`1e-12`) o relativo a `|field|²` local?
- [ ] ¿Los tests de edge cases que setean `bundle.sx/sy` deberían eliminarse o reescribirse para ejercitar slopes reales?