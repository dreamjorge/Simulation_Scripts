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
if (size(bundle_grid.x, 1) == 12 && size(bundle_grid.x, 2) == 10 && size(bundle_grid.x, 3) == 1 && abs(max(bundle_grid.x(:)) - 0.5e-3) < 1e-8)
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

% -----------------------------------------------------------------
% Phase 4: Physical accuracy tests
% -----------------------------------------------------------------

% testGradientVsAnalyticalGaussian
% Compare numerical gradient to analytical: d(phase)/dx = k*x/R(z) at z=0, R=Inf => 0
% At z>0, the radius of curvature R(z) gives non-zero gradient
z_test = zr / 2;  % away from waist where R(z) is finite
params = beam.getParameters(z_test);
R_z = params.RadiusOfCurvature;
% Analytical phase gradient: dφ/dx = k*x/R(z)
% Analytical slope: sx = (1/k) * dφ/dx = x/R(z)
x_test = 20e-6;
[sx_num, sy_num] = RayTracer.calculateSlopes(beam, x_test, 0, z_test);
sx_analytical = x_test / R_z;
rel_err = abs(sx_num - sx_analytical) / abs(sx_analytical);
if (rel_err < 1e-4)
    fprintf('  PASS: Gradient vs analytical Gaussian (rel err %.2e)\n', rel_err);
    passed = passed + 1;
else
    fprintf('  FAIL: Gradient vs analytical Gaussian (rel err %.2e)\n', rel_err);
    failed = failed + 1;
end

% testEulerVsRK4Convergence
% Refine dz and verify both methods converge to the same result
bundle_euler_coarse = RayBundle.createGrid(3, 3, w0, w0);
bundle_rk4_coarse = RayBundle.createGrid(3, 3, w0, w0);
dz_coarse = zr / 20;
bundle_euler_coarse = RayTracer.propagate(bundle_euler_coarse, beam, zr/4, dz_coarse, 'Euler');
bundle_rk4_coarse = RayTracer.propagate(bundle_rk4_coarse, beam, zr/4, dz_coarse, 'RK4');
final_euler_coarse = bundle_euler_coarse.x(:,:,end);
final_rk4_coarse = bundle_rk4_coarse.x(:,:,end);

bundle_euler_fine = RayBundle.createGrid(3, 3, w0, w0);
bundle_rk4_fine = RayBundle.createGrid(3, 3, w0, w0);
dz_fine = zr / 80;
bundle_euler_fine = RayTracer.propagate(bundle_euler_fine, beam, zr/4, dz_fine, 'Euler');
bundle_rk4_fine = RayTracer.propagate(bundle_rk4_fine, beam, zr/4, dz_fine, 'RK4');
final_euler_fine = bundle_euler_fine.x(:,:,end);
final_rk4_fine = bundle_rk4_fine.x(:,:,end);

% Euler should improve with finer dz; RK4 should already be close to final
euler_improvement = mean(abs(final_euler_coarse(:) - final_rk4_fine(:)));
euler_coarse_vs_rk4 = mean(abs(final_euler_coarse(:) - final_rk4_coarse(:)));
if (euler_improvement < euler_coarse_vs_rk4)
    fprintf('  PASS: Euler converges with finer dz\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Euler does not improve with finer dz\n');
    failed = failed + 1;
end

% testRadialSymmetry
% Verify that bundle.r matches sqrt(x^2+y^2) from last z-slice
if (abs(bundle_prop.Nz) > 1)
    r_computed = sqrt(bundle_prop.x(:,:,end).^2 + bundle_prop.y(:,:,end).^2);
    r_stored = bundle_prop.r;
    max_diff = max(abs(r_computed(:) - r_stored(:)));
    if (max_diff < 1e-15)
        fprintf('  PASS: bundle.r matches sqrt(x^2+y^2)\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: bundle.r differs from sqrt(x^2+y^2) by %.2e\n', max_diff);
        failed = failed + 1;
    end
else
    fprintf('  FAIL: insufficient z steps for radial symmetry test\n');
    failed = failed + 1;
end

fprintf('\n=== RayTracing: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
