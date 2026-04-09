# Simulation_Scripts

Scripts for simulation of optical beam propagation (Gaussian, Hermite-Gauss, Laguerre-Gauss beams).

Author: Ugalde-Ontiveros J.A.

## Estructura

```
Simulation_Scripts/
├── ParaxialBeams/
│   ├── GaussianBeam.m
│   ├── HermiteBeam.m
│   ├── LaguerreBeam.m
│   ├── ElegantHermiteBeam.m
│   ├── ElegantLaguerreBeam.m
│   ├── GaussianParameters.m
│   ├── HermiteParameters.m
│   ├── LaguerreParameters.m
│   ├── PhysicalConstants.m
│   ├── GridUtils.m
│   ├── FFTUtils.m
│   ├── RayBundle.m / RayTracer.m
│   └── Addons/
├── examples/
│   ├── MainGauss_refactored.m
│   └── MainMultiMode.m
├── tests/
│   ├── test_all.m
│   └── portable_runner.m
└── ExampleRayTracing.m
```

## ⭐ Nuevas Clases (v2.0)

### PhysicalConstants
Constantes físicas y métodos utilitarios:
```matlab
PC = PhysicalConstants;
k = PC.waveNumber(lambda);
zr = PC.rayleighDistance(w0, lambda);
w = PC.waistAtZ(w0, z, lambda, zr);
R = PC.radiusOfCurvature(z, zr);
gouy = PC.gouyPhase(z, zr);
```

### GridUtils
Generación de grids computacionales:
```matlab
grid = GridUtils(Nx, Ny, Dx, Dy);
[X, Y] = grid.create2DGrid();
[Kx, Ky] = grid.createFreqGrid();
[r, theta] = grid.createPolarGrid();
```

### FFTUtils
Operaciones FFT normalizadas:
```matlab
fftOps = FFTUtils(true, true);  % normalize, shift
G = fftOps.fft2(g);              % FFT hacia adelante
g = fftOps.ifft2(G);              % FFT inversa
H = fftOps.transferFunction(kx, ky, z, lambda);  % Función de transferencia
gProp = fftOps.propagate(g, kx, ky, z, lambda);  % Propagación angular
```

## Uso

## Beam API Contract

- canonical field entrypoint: `opticalField(...)`
- `Parameters` define constantes/modelo, no reemplaza argumentos dinámicos ambiguamente
- cada beam debe declarar sus coordenadas aceptadas
- cualquier desvío temporal queda documentado hasta el post-merge cleanup

### Canonical Examples

- `examples/MainGauss_refactored.m` - entrada principal para propagación Gaussiana con utilidades modernas
- `examples/MainMultiMode.m` - entrada recomendada para modos Hermite y Laguerre
- `ExampleRayTracing.m` - entrada recomendada para ray tracing integrado

### MATLAB
```matlab
addpath ParaxialBeams
addpath ParaxialBeams\Addons

% Usar constantes físicas
PC = PhysicalConstants;
zr = PC.rayleighDistance(100e-6, 632.8e-9);

% Crear grid
grid = GridUtils(1024, 1024, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

% FFT normalizado
fftOps = FFTUtils;
G = fftOps.fft2(field);
```

### Compatibilidad

- Octave: verificado localmente con `tests/test_all.m`
- MATLAB completo: requerido para validar `matlab -batch "run('tests/test_all.m')"`
- MATLAB Runtime solo no alcanza para ejecutar la suite del repo

## Referencias

- Kogelnik, H., & Li, T. (1966). Laser beams and resonators. Applied Optics.
- Siegman, A. E. (1986). Lasers. University Science Books.
