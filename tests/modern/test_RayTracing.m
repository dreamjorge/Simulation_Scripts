% Compatible with GNU Octave and MATLAB
% Tests for Ray tracing integration (script-style for portable runner)

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

fprintf('=== RayTracing Tests ===\n\n');
passed = 0;
failed = 0;

% -----------------------------------------------------------------
% RayBundle initialization
% -----------------------------------------------------------------

% testRayBundleGridInitialization
bundle_grid = RayBundle.createGrid(10, 12, 1e-3, 1.2e-3);
if (isequal(size(bundle_grid.x), [12, 10, 1]) && abs(max(bundle_grid.x(:)) - 0.5e-3) < 1e-8)
    fprintf('  PASS: RayBundle createGrid initialization\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle createGrid initialization\n');
    failed = failed + 1;
end

% testRayBundleConcentricInitialization
bundle_conc = RayBundle.createConcentric(5, 8, 1e-3);
if (bundle_conc.Nx == 5 && bundle_conc.Ny == 8)
    fprintf('  PASS: RayBundle createConcentric initialization\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle createConcentric initialization\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% RayTracer slope calculations
% -----------------------------------------------------------------

lambda = 632.8e-9;
w0 = 100e-6;
beam = GaussianBeam(w0, lambda);

% testRayTracerSlopesAtCenter
[sx0, sy0] = RayTracer.calculateSlopes(beam, 0, 0, 0);
if (abs(sx0) < 1e-10 && abs(sy0) < 1e-10)
    fprintf('  PASS: RayTracer slopes at center\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracer slopes at center\n');
    failed = failed + 1;
end

% testRayTracerSlopesOffCenterAtZ0
[sx1, sy1] = RayTracer.calculateSlopes(beam, 50e-6, 0, 0);
if (abs(sx1) < 1e-6 && abs(sy1) < 1e-6)
    fprintf('  PASS: RayTracer slopes off-center at z=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracer slopes off-center at z=0\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% Free-space propagation
% -----------------------------------------------------------------

zr = pi * w0^2 / lambda;
bundle_prop = RayBundle.createGrid(5, 5, 200e-6, 200e-6);
dz = zr / 10;
bundle_prop = RayTracer.propagate(bundle_prop, beam, zr, dz, 'RK4');

% testPropagationAddsSteps
if (bundle_prop.Nz > 1)
    fprintf('  PASS: RayTracer propagation adds z steps\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracer propagation adds z steps\n');
    failed = failed + 1;
end

% testPropagationDivergence
last_x = bundle_prop.x(1, end, end);
first_x = bundle_prop.x(1, end, 1);
if (abs(last_x) > abs(first_x))
    fprintf('  PASS: RayTracer propagation diverges rays\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracer propagation diverges rays\n');
    failed = failed + 1;
end

fprintf('\n=== RayTracing: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
