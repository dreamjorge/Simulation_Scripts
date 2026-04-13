# Refactor Completo — `refactor/utility-classes` Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Arreglar todos los issues técnicos del branch refactor/utility-classes antes de hacer PR.

**Architecture:**
- **GridUtils**: Fix static methods `meshgrid2D` y `freqGrid` para soportar Nx/Ny separados
- **ElegantHermiteBeam**: Eliminar uso innecesario de `cart2pol` — compute R directamente
- **Elegant beam naming**: Estandarizar a lowercase `x,y` para coordenadas cartesianas
- **Backwards-compat wrappers**: Remover `HermiteBeam.hermitePoly()` y `LaguerreParameters.getAssociatedLaguerrePolynomial()`
- **HankelLaguerre**: Implementar con formula real descubierta en git history
- **AnalysisUtils.combinedHankelWave**: Decidir si implementar o remover stub
- **Addons/Addons que llaman a stubs**: getPropagateCylindricalRays y getPropagateRay — SI se usan
- **Test compatibility**: Reemplazar `exit(0)/exit(1)` con `assert()` para MATLAB compatibilidad

**Tech Stack:** Octave/MATLAB, classdef, (no framework — vanilla scripts)

---

## Task 1: Investigar HankelLaguerre y archivos relacionados ✅ (COMPLETADO)

**Investigación completada:**

### Archivos que usan getPropagateCylindricalRays:
- `MainLaguerre.m` — líneas 206, 218
- `MainLaguerre2.m` — líneas 203, 215
- `MaineLaguerreForTest.m` — líneas 203, 215
- `AnalysisExperimentalLaguerre.m` — líneas 207-225

### Fórmula real de HankelLaguerre (descubierta en commit 981b7b1):

```matlab
% Del @HankelLaguerre legacy (MATLAB class folder):
if nh == 1
    Hankel.OpticalFieldLaguerre = LB.OpticalFieldLaguerre + 1i*XLB.OpticalFieldLaguerre;
elseif nh == 2
    Hankel.OpticalFieldLaguerre = LB.OpticalFieldLaguerre - 1i*XLB.OpticalFieldLaguerre;
end
```

Donde XLaguerreBeam usa:
- `exp(-1i*abs(p)*theta)` vs `exp(1i*l*theta)` en fase
- `XAssociatedLaguerrePolynomial(l,abs(p),xArg)` = `associatedLaguerre(p, l, x)` (orden swapped!)

---

## Task 2: Implementar HankelLaguerre — Formula real

**Files:**
- Modify: `ParaxialBeams/HankelLaguerre.m`
- Create: `tests/test_HankelLaguerre.m`

**Step 1: Escribir test primero**

```matlab
#!/usr/bin/env octave
% Tests for HankelLaguerre

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== HankelLaguerre Tests ===\n\n');
passed = 0;
failed = 0;

N = 64; D = 1e-3;
grid = GridUtils(N, N, D, D);
[X, Y] = grid.create2DGrid();
[r, theta] = cart2pol(X, Y);

params = LaguerreParameters(0, 1e-3, 632e-9, 1, 0);

% testConstructorHankelType1
try
    HL = HankelLaguerre(r, theta, params, 1);
    if ~isempty(HL.OpticalFieldLaguerre) && size(HL.OpticalFieldLaguerre) == size(r)
        fprintf('  PASS: constructor hankelType=1\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: OpticalFieldLaguerre empty or wrong size\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testConstructorHankelType2
try
    HL = HankelLaguerre(r, theta, params, 2);
    if ~isempty(HL.OpticalFieldLaguerre)
        fprintf('  PASS: constructor hankelType=2\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: OpticalFieldLaguerre empty\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testConstructorDefaultType
try
    HL = HankelLaguerre(r, theta, params);
    if ~isempty(HL.OpticalFieldLaguerre)
        fprintf('  PASS: constructor default type\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default type failed\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testHankelType1PlusType2
try
    HL1 = HankelLaguerre(r, theta, params, 1);
    HL2 = HankelLaguerre(r, theta, params, 2);
    % H1 + H2 should equal 2*LB (imag parts cancel)
    diff = HL1.OpticalFieldLaguerre + HL2.OpticalFieldLaguerre;
    LB = LaguerreBeam(r, theta, params);
    if max(max(abs(diff - 2*LB.OpticalField))) < 1e-10
        fprintf('  PASS: H1 + H2 = 2*LB\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: H1 + H2 != 2*LB\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

fprintf('\n=== HankelLaguerre: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
```

Run: `octave tests/test_HankelLaguerre.m`
Expected: FAIL (current stub throws error)

**Step 2: Implement HankelLaguerre**

Reemplazar todo el contenido de `ParaxialBeams/HankelLaguerre.m`:

```matlab
classdef HankelLaguerre
    % HankelLaguerre - Hankel-type Laguerre-Gaussian beam field
    % Used in Hankel-based ray tracing (HankelLaguerrePropagation.m, etc.)
    %
    % Implements: HL^{(p,l)}_1 = LB + 1i*XLG
    %             HL^{(p,l)}_2 = LB - 1i*XLG
    %
    % Where LB is standard LaguerreBeam and XLG is the Hilbert-transformed
    % counterpart with:
    %   - Phase: exp(-1i*p*theta) instead of exp(1i*l*theta)
    %   - Polynomial: AssociatedLaguerre(p, l, x) (args in different order)

    properties
        OpticalFieldLaguerre  % Complex field array
    end

    methods
        function obj = HankelLaguerre(r, theta, params, hankelType)
            % Constructor
            % r: radial coordinate matrix
            % theta: angular coordinate matrix
            % params: LaguerreParameters object
            % hankelType: 1 (H^(1)) or 2 (H^(2))
            
            if nargin < 4
                hankelType = 1;
            end
            
            obj.OpticalFieldLaguerre = computeHankelField(r, theta, params, hankelType);
        end
    end
end

function field = computeHankelField(r, theta, params, hankelType)
    l = params.l;
    p = params.p;
    w = params.Waist;
    
    % Standard LG amplitude term: (sqrt(2)*r/w)^|l|
    amp = (sqrt(2) * r ./ w).^abs(l);
    
    % x argument for Laguerre polynomial
    xArg = 2 * r.^2 ./ w.^2;
    
    % Standard LG polynomial L_p^l(x)
    Lpl = PolynomialUtils.associatedLaguerre(p, l, xArg);
    
    % Gaussian carrier field
    GB = GaussianBeam(r, params);
    GField = GB.OpticalField;
    
    % Standard LG field (LB)
    LB_field = amp .* Lpl .* exp(1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;
    
    % XLG (Hilbert-transformed) field:
    % - Phase uses exp(-1i*p*theta) instead of exp(1i*l*theta)
    % - Polynomial uses associatedLaguerre(p, l, x) which is what we already have!
    XLG_field = amp .* Lpl .* exp(-1i * p * theta) .* exp(1i * params.PhiPhase) .* GField;
    
    % Combine via Hankel type
    if hankelType == 1
        field = LB_field + 1i * XLG_field;
    else
        field = LB_field - 1i * XLG_field;
    end
end
```

Run test: `octave tests/test_HankelLaguerre.m`
Expected: PASS

**Step 3: Commit**
```bash
git add ParaxialBeams/HankelLaguerre.m tests/test_HankelLaguerre.m
git commit -m "feat(HankelLaguerre): implement Hankel-type LG beam field"
```

---

## Task 3: Fix GridUtils static methods

**Files:**
- Modify: `ParaxialBeams/GridUtils.m`
- Modify: `tests/test_GridUtils.m`

**Step 1: Fix meshgrid2D static method**

```matlab
% ANTES (line 72-75):
function [X, Y] = meshgrid2D(N, D)
    n = -N/2:N/2-1;
    x = n * (D / N);
    [X, Y] = meshgrid(x, x);
end

% NUEVO (line 72-79):
function [X, Y] = meshgrid2D(Nx, Ny, Dx, Dy)
    % meshgrid2D - Create 2D coordinate grid
    % Supports asymmetric grids when Nx~=Ny or Dx~=Dy
    if nargin < 3
        % Backwards compat: single N, D expands to Nx=Ny, Dx=Dy
        Dy = Ny; Dx = Ny; Ny = Nx;
    end
    nx = -Nx/2:Nx/2-1;
    ny = -Ny/2:Ny/2-1;
    x = nx * (Dx / Nx);
    y = ny * (Dy / Ny);
    [X, Y] = meshgrid(x, y);
end
```

**Step 2: Fix freqGrid static method**

```matlab
% ANTES (line 81-85):
function [Kx, Ky] = freqGrid(N, D)
    n = -N/2:N/2-1;
    u = n * (1 / D);
    [Kx, Ky] = meshgrid(2*pi*u, 2*pi*u);
end

% NUEVO (line 81-91):
function [Kx, Ky] = freqGrid(Nx, Ny, Dx, Dy)
    % freqGrid - Create frequency domain grid
    % Supports asymmetric grids when Nx~=Ny or Dx~=Dy
    if nargin < 4
        % Backwards compat: single N, D expands to Nx=Ny, Dx=Dy
        Dy = Dx; Ny = Nx;
    end
    nx = -Nx/2:Nx/2-1;
    ny = -Ny/2:Ny/2-1;
    u = nx * (1 / Dx);
    v = ny * (1 / Dy);
    [Kx, Ky] = meshgrid(2*pi*u, 2*pi*v);
end
```

**Step 3: Fix polarGrid (llama a meshgrid2D)**

```matlab
% ANTES (line 93-96):
function [r, theta] = polarGrid(N, D)
    [X, Y] = GridUtils.meshgrid2D(N, D);
    [r, theta] = cart2pol(X, Y);
end

% NUEVO (line 93-101):
function [r, theta] = polarGrid(Nx, Ny, Dx, Dy)
    % polarGrid - Create polar coordinate grid
    % Supports asymmetric grids when Nx~=Ny or Dx~=Dy
    if nargin < 4
        % Backwards compat: single N, D expands to Nx=Ny, Dx=Dy
        Dy = Dx; Ny = Nx;
    end
    [X, Y] = GridUtils.meshgrid2D(Nx, Ny, Dx, Dy);
    [r, theta] = cart2pol(X, Y);
end
```

**Step 4: Correr tests**

Run: `octave tests/test_GridUtils.m`
Expected: PASS (backwards compat mantener N,D funcionando)

**Step 5: Commit**
```bash
git add ParaxialBeams/GridUtils.m tests/test_GridUtils.m
git commit -m "fix(GridUtils): static methods now support asymmetric grids"
```

---

## Task 4: Estandarizar naming en Elegant beams

**Files:**
- Modify: `ParaxialBeams/ElegantHermiteBeam.m`

**Step 1: Analizar convencion actual**

- `HermiteBeam`: usa `x`, `y` (lowercase)
- `ElegantHermiteBeam`: usa `X`, `Y` (uppercase) — INCONSISTENTE
- `ElegantLaguerreBeam`: usa `r`, `theta` (correcto para cilíndricas)

Elegir convencion: **lowercase** (`x`, `y`) para consistencia con `HermiteBeam`.

**Step 2: Fix ElegantHermiteBeam**

```matlab
% ANTES (propiedades, lines 7-8):
X               % X coordinate matrix
Y               % Y coordinate matrix

% NUEVO (propiedades):
x               % x coordinate matrix
y               % y coordinate matrix

% Constructor (lines 12-17) cambiar X,Y por x,y en todo
% Tambien line 21: argX = sqrt_alpha .* X; -> argX = sqrt_alpha .* x;
% Line 22: argY = sqrt_alpha .* Y; -> argY = sqrt_alpha .* y;
% Line 29: [R, ~] = cart2pol(X, Y); -> [R, ~] = cart2pol(x, y);
```

**Step 3: Commit**
```bash
git add ParaxialBeams/ElegantHermiteBeam.m
git commit -m "refactor(ElegantHermiteBeam): standardize to lowercase x,y"
```

---

## Task 5: Eliminar backwards-compat wrappers

**Files:**
- Modify: `ParaxialBeams/HermiteBeam.m` (remover método estático hermitePoly)
- Modify: `ParaxialBeams/LaguerreParameters.m` (remover método estático getAssociatedLaguerrePolynomial)

**Step 1: Buscar usages de estos wrappers**

```bash
grep -r "HermiteBeam.hermitePoly\|LaguerreParameters.getAssociatedLaguerrePolynomial" /root/Simulation_Scripts --include="*.m"
```

Expected: solo las definiciones (no se usan)

**Step 2: Remover wrappers**

De `HermiteBeam.m` líneas 45-50 remover el bloque completo:
```matlab
methods (Static)
    function H = hermitePoly(n, x)
        H = PolynomialUtils.hermitePoly(n, x);
    end
end
```

De `LaguerreParameters.m` líneas 59-62 remover solo `getAssociatedLaguerrePolynomial`:
```matlab
methods (Static)
    function wL = getWaist(z, w0, zr, l, p)
        w = w0 * sqrt(1 + (z/zr).^2);
        wL = w * sqrt(2*p + abs(l) + 1);
    end
    
    function L = getAssociatedLaguerrePolynomial(p, l, x)
        L = PolynomialUtils.associatedLaguerre(p, l, x);
    end
end
```

**Mantener** `getWaist` (es útil) pero remover `getAssociatedLaguerrePolynomial`.

**Step 3: Commit**
```bash
git add ParaxialBeams/HermiteBeam.m ParaxialBeams/LaguerreParameters.m
git commit -m "refactor: remove backwards-compat wrappers"
```

---

## Task 6: Fix ElegantHermiteBeam cart2pol innecesario

**Files:**
- Modify: `ParaxialBeams/ElegantHermiteBeam.m`

**Step 1: Analizar**

El código actual (línea 29):
```matlab
[R, ~] = cart2pol(X, Y);
GB = GaussianBeam(R, params);
```

Pero `GaussianBeam` recibe `r` (radial), no necesita `theta`. El `cart2pol` computa theta que se descarta.

**Step 2: Compute R directamente**

```matlab
% ANTES:
[R, ~] = cart2pol(X, Y);
GB = GaussianBeam(R, params);

% NUEVO:
R = sqrt(X.^2 + Y.^2);
GB = GaussianBeam(R, params);
```

**Step 3: Commit**
```bash
git add ParaxialBeams/ElegantHermiteBeam.m
git commit -m "fix(ElegantHermiteBeam): compute radius directly instead of cart2pol"
```

---

## Task 7: Test compatibility — reemplazar exit() con assert()

**Files:**
- Modify: `tests/*.m` (todos los archivos de test)

**Step 1: Reemplazar exit() con assert()**

Patrón a buscar en cada test file:
```matlab
% ANTES:
if failed == 0
    exit(0);
else
    exit(1);
end

% NUEVO:
if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
```

Hacerlo para todos los test files:
- test_all.m
- test_GridUtils.m
- test_FFTUtils.m
- test_GaussianParameters.m
- test_PhysicalConstants.m
- test_HermiteBeam.m
- test_AnalysisUtils.m
- test_ElegantLaguerreBeam.m
- test_ElegantHermiteBeam.m
- test_LaguerreBeam.m
- test_GaussianBeam.m
- test_ElegantLaguerreParameters.m
- test_ElegantHermiteParameters.m
- test_HermiteParameters.m
- test_LaguerreParameters.m

**Step 2: Commit**
```bash
git add tests/*.m
git commit -m "test: replace exit() with assert() for MATLAB compatibility"
```

---

## Task 8: Revisar AnalysisUtils.combinedHankelWave

**Files:**
- Analyze: `ParaxialBeams/AnalysisUtils.m`
- Analyze: `ParaxialBeams/HankelLaguerre.m` (note sobre combinedHankelWave)

**Step 1: Buscar usages**

```bash
grep -r "combinedHankelWave" /root/Simulation_Scripts --include="*.m"
```

**Step 2: Si no se usa**

Este stub hace referencia a `HankelHermiteSlices.m`. No se encontró en el código actual. Remover el stub y su comentario de error.

---

## Task 9: Run full test suite y verification

**Step 1: Correr todos los tests**

```bash
cd /root/Simulation_Scripts && octave tests/test_all.m
```

Expected: Todos PASS

**Step 2: Verificar coverage**

```bash
octave tests/run_coverage.m
```

---

## Pre-PR Checklist

- [ ] HankelLaguerre implementado y testeado
- [ ] GridUtils static methods fix con backwards compat
- [ ] ElegantHermiteBeam naming (x,y) y cart2pol fix
- [ ] Backwards-compat wrappers removidos
- [ ] Tests con error() en vez de exit()
- [ ] AnalysisUtils.combinedHankelWave stub removido si no se usa
- [ ] Todos los tests pasando
- [ ] No nuevas advertencias

---

## Riscos

1. **Implementación HankelLaguerre**: La formula real requiere XLaguerreBeam que tiene normalización extra — verificar con scripts reales
2. **Breaking backwards compat**: GridUtils ahora soporta Nx,Ny,Dx,Dy — signature changed pero con backwards compat
3. **Breaking backwards compat wrappers**: Los wrappers removidos eran solo para backwards — verificar que no se usen externamente

**Mitigación**: Tests primero, commits atômicos, validación manual antes de PR.