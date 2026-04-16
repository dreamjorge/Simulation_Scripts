% Compatible with GNU Octave and MATLAB
% Tests for HankelRayTracer and HankelRayTracePropagator

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== HankelRayTracing Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
zr = pi * w0^2 / lambda;

% -----------------------------------------------------------------
% RayBundle hankelType tracking
% -----------------------------------------------------------------

% testRayBundleHankelTypeDefault
bundle = RayBundle.createGrid(4, 4, 1e-4, 1e-4);
if all(bundle.ht(:) == 1)
    fprintf('  PASS: RayBundle ht defaults to 1\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle ht default\n');
    failed = failed + 1;
end

% testRayBundleHankelTypeSet
bundle.ht(:) = 11;
if all(bundle.ht(:) == 11)
    fprintf('  PASS: RayBundle ht assignment\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle ht assignment\n');
    failed = failed + 1;
end

% testRayBundleAddStepPreservesHt
bundle.addStep(bundle.x(:,:,1), bundle.y(:,:,1), bundle.z(:,:,1) + 1e-5, zeros(4), zeros(4));
if all(bundle.ht(:,:,end) == 11)
    fprintf('  PASS: addStep preserves ht\n');
    passed = passed + 1;
else
    fprintf('  FAIL: addStep preserves ht\n');
    failed = failed + 1;
end

% testRayBundleAddStepExplicitHt
bundle.addStep(bundle.x(:,:,1), bundle.y(:,:,1), bundle.z(:,:,1) + 2e-5, zeros(4), zeros(4), 22 * ones(4));
if all(bundle.ht(:,:,end) == 22)
    fprintf('  PASS: addStep explicit ht\n');
    passed = passed + 1;
else
    fprintf('  FAIL: addStep explicit ht\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% HankelRayTracer with HankelHermite
% -----------------------------------------------------------------

% testHankelRayTracerHermite
beam_h = HankelHermite(w0, lambda, 1, 1, 11);
bundle_h = RayBundle.createGrid(3, 3, w0, w0);
bundle_h.ht(:) = 11;
try
    bundle_h = HankelRayTracer.propagate(bundle_h, beam_h, zr/10, zr/100, 'RK4');
    nSteps = size(bundle_h.x, 3);
    if nSteps > 1 && all(all(isfinite(bundle_h.x(:,:,end)))) && all(all(isfinite(bundle_h.y(:,:,end))))
        fprintf('  PASS: HankelRayTracer Hermite propagation\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: HankelRayTracer Hermite propagation (infinite coords)\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: HankelRayTracer Hermite (%s)\n', ME.message);
    failed = failed + 1;
end

% testHankelRayTracerHermiteDiverges
if nSteps > 1
    x_start = bundle_h.x(:,:,1);
    x_end = bundle_h.x(:,:,end);
    if max(abs(x_end(:))) >= max(abs(x_start(:)))
        fprintf('  PASS: Hermite rays diverge\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: Hermite rays did not diverge\n');
        failed = failed + 1;
    end
else
    fprintf('  FAIL: Hermite rays no steps\n');
    failed = failed + 1;
end

% testHankelRayTracerHermiteHtPreserved
if all(bundle_h.ht(:,:,end) == 11)
    fprintf('  PASS: Hermite ht preserved (no axis crossing)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Hermite ht changed\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% HankelRayTracer with HankelLaguerre
% -----------------------------------------------------------------

% testHankelRayTracerLaguerre
beam_l = HankelLaguerre(w0, lambda, 1, 0, 1);
bundle_l = RayBundle.createConcentric(3, 6, w0);
bundle_l.ht(:) = 1;
try
    bundle_l = HankelRayTracer.propagate(bundle_l, beam_l, zr/10, zr/100, 'RK4');
    nSteps_l = size(bundle_l.x, 3);
    if nSteps_l > 1 && all(all(isfinite(bundle_l.x(:,:,end))))
        fprintf('  PASS: HankelRayTracer Laguerre propagation\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: HankelRayTracer Laguerre propagation\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: HankelRayTracer Laguerre (%s)\n', ME.message);
    failed = failed + 1;
end

% -----------------------------------------------------------------
% HankelRayTracePropagator (Strategy pattern)
% -----------------------------------------------------------------

% testHankelRayTracePropagatorHermite
grid = GridUtils(4, 4, 2*w0, 2*w0);
prop_h = HankelRayTracePropagator(grid, 11, 'RK4', zr/50);
try
    bundle_ph = prop_h.propagate(beam_h, zr/10);
    if size(bundle_ph.x, 3) > 1 && all(bundle_ph.ht(:,:,1) == 11)
        fprintf('  PASS: HankelRayTracePropagator Hermite\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: HankelRayTracePropagator Hermite\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: HankelRayTracePropagator Hermite (%s)\n', ME.message);
    failed = failed + 1;
end

% testHankelRayTracePropagatorLaguerre
prop_l = HankelRayTracePropagator(grid, 1, 'RK4', zr/50);
try
    bundle_pl = prop_l.propagate(beam_l, zr/10);
    if size(bundle_pl.x, 3) > 1
        fprintf('  PASS: HankelRayTracePropagator Laguerre\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: HankelRayTracePropagator Laguerre\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: HankelRayTracePropagator Laguerre (%s)\n', ME.message);
    failed = failed + 1;
end

% testEulerMethod
prop_euler = HankelRayTracePropagator(grid, 11, 'Euler', zr/50);
try
    bundle_euler = prop_euler.propagate(beam_h, zr/10);
    if size(bundle_euler.x, 3) > 1
        fprintf('  PASS: Euler method works\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: Euler method\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: Euler method (%s)\n', ME.message);
    failed = failed + 1;
end

fprintf('\n=== HankelRayTracing: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
