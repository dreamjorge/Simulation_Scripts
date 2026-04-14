% Compatible with GNU Octave and MATLAB
% Tests for IPropagator, AnalyticPropagator, FFTPropagator, RayTracePropagator

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== Propagators Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
beam = GaussianBeam(w0, lambda);

% -----------------------------------------------------------------
% AnalyticPropagator
% -----------------------------------------------------------------

% testAnalyticConstructor
ap = AnalyticPropagator(grid);
if (isa(ap, 'IPropagator') && isa(ap, 'AnalyticPropagator'))
    fprintf('  PASS: AnalyticPropagator constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator constructor\n');
    failed = failed + 1;
end

% testAnalyticPropagateSize
field = ap.propagate(beam, 0.1);
if (size(field) == [64, 64])
    fprintf('  PASS: AnalyticPropagator field size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator field size\n');
    failed = failed + 1;
end

% testAnalyticPropagateFinite
if (all(all(isfinite(field))))
    fprintf('  PASS: AnalyticPropagator finite output\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator finite output\n');
    failed = failed + 1;
end

% testAnalyticMatchesDirectCall
[X, Y] = grid.create2DGrid();
field_direct = beam.opticalField(X, Y, 0.1);
if (max(max(abs(field - field_direct))) < 1e-12)
    fprintf('  PASS: AnalyticPropagator matches direct opticalField\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator vs direct opticalField\n');
    failed = failed + 1;
end

% testAnalyticAtZ0
field_z0 = ap.propagate(beam, 0);
if (all(all(isfinite(field_z0))))
    fprintf('  PASS: AnalyticPropagator at z=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator at z=0\n');
    failed = failed + 1;
end

% testAnalyticWithHermiteBeam
hb = HermiteBeam(w0, lambda, 1, 1);
field_h = ap.propagate(hb, 0.05);
if (all(all(isfinite(field_h))) && size(field_h) == [64, 64])
    fprintf('  PASS: AnalyticPropagator with HermiteBeam\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator with HermiteBeam\n');
    failed = failed + 1;
end

% testAnalyticWithLaguerreBeam
lb = LaguerreBeam(w0, lambda, 1, 0);
field_l = ap.propagate(lb, 0.05);
if (all(all(isfinite(field_l))) && size(field_l) == [64, 64])
    fprintf('  PASS: AnalyticPropagator with LaguerreBeam\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator with LaguerreBeam\n');
    failed = failed + 1;
end

% testAnalyticZ0Property
ap_z1 = AnalyticPropagator(grid, 0.01);
if (ap_z1.z0 == 0.01)
    fprintf('  PASS: AnalyticPropagator z0 property\n');
    passed = passed + 1;
else
    fprintf('  FAIL: AnalyticPropagator z0 property\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% FFTPropagator
% -----------------------------------------------------------------

% testFFTConstructor
fp = FFTPropagator(grid, lambda);
if (isa(fp, 'IPropagator') && isa(fp, 'FFTPropagator'))
    fprintf('  PASS: FFTPropagator constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator constructor\n');
    failed = failed + 1;
end

% testFFTPropagateSize
field_fft = fp.propagate(beam, 0.1);
if (size(field_fft) == [64, 64])
    fprintf('  PASS: FFTPropagator field size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator field size\n');
    failed = failed + 1;
end

% testFFTPropagateFinite
if (all(all(isfinite(field_fft))))
    fprintf('  PASS: FFTPropagator finite output\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator finite output\n');
    failed = failed + 1;
end

% testFFTvsAnalyticRelativeError
% For a Gaussian beam in the paraxial regime, FFT and analytic should agree
% closely at center (relative error < 5% is acceptable for 64x64 grid)
field_analytic = ap.propagate(beam, 0.1);
amp_fft     = abs(field_fft(33, 33));
amp_analytic = abs(field_analytic(33, 33));
rel_err = abs(amp_fft - amp_analytic) / amp_analytic;
if (rel_err < 0.05)
    fprintf('  PASS: FFT vs analytic < 5%% at center\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFT vs analytic rel_err=%.4f\n', rel_err);
    failed = failed + 1;
end

% testFFTAtZ0Identity
field_fft_z0 = fp.propagate(beam, 0);
field_ref_z0 = ap.propagate(beam, 0);
rel_err_z0 = max(max(abs(abs(field_fft_z0) - abs(field_ref_z0)))) / max(max(abs(field_ref_z0)));
if (rel_err_z0 < 0.01)
    fprintf('  PASS: FFTPropagator at z=0 near identity\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator at z=0 rel_err=%.4f\n', rel_err_z0);
    failed = failed + 1;
end

% testFFTLambdaStored
if (fp.Lambda == lambda)
    fprintf('  PASS: FFTPropagator Lambda stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator Lambda stored\n');
    failed = failed + 1;
end

% testFFTz0Property
fp_z1 = FFTPropagator(grid, lambda, 0.01);
if (fp_z1.z0 == 0.01)
    fprintf('  PASS: FFTPropagator z0 property\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator z0 property\n');
    failed = failed + 1;
end

% testFFTWithHermiteBeam
field_fft_h = fp.propagate(HermiteBeam(w0, lambda, 1, 1), 0.05);
if (all(all(isfinite(field_fft_h))) && size(field_fft_h) == [64, 64])
    fprintf('  PASS: FFTPropagator with HermiteBeam\n');
    passed = passed + 1;
else
    fprintf('  FAIL: FFTPropagator with HermiteBeam\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% RayTracePropagator
% -----------------------------------------------------------------

% testRayTraceConstructor
rp = RayTracePropagator(grid);
if (isa(rp, 'IPropagator') && isa(rp, 'RayTracePropagator'))
    fprintf('  PASS: RayTracePropagator constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator constructor\n');
    failed = failed + 1;
end

% testRayTraceDefaultMethod
if (strcmp(rp.Method, 'RK4'))
    fprintf('  PASS: RayTracePropagator default method RK4\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator default method\n');
    failed = failed + 1;
end

% testRayTraceReturnsBundle
% Use a coarse grid to keep it fast
grid_small = GridUtils(8, 8, 1e-3, 1e-3);
rp_small   = RayTracePropagator(grid_small, 'Euler', 0.005);
bundle = rp_small.propagate(beam, 0.01);
if (isa(bundle, 'RayBundle'))
    fprintf('  PASS: RayTracePropagator returns RayBundle\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator returns RayBundle\n');
    failed = failed + 1;
end

% testRayTraceBundleHasSteps
if (bundle.Nz > 1)
    fprintf('  PASS: RayBundle has multiple z steps\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle Nz <= 1\n');
    failed = failed + 1;
end

% testRayTraceBundleFinite
if (all(all(all(isfinite(bundle.x)))) && all(all(all(isfinite(bundle.y)))))
    fprintf('  PASS: RayBundle coordinates finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle coordinates not finite\n');
    failed = failed + 1;
end

% testRayTraceMethodStored
rp_euler = RayTracePropagator(grid, 'Euler');
if (strcmp(rp_euler.Method, 'Euler'))
    fprintf('  PASS: RayTracePropagator method stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator method stored\n');
    failed = failed + 1;
end

% testRayTraceDzStored
rp_dz = RayTracePropagator(grid, 'RK4', 1e-3);
if (rp_dz.dz == 1e-3)
    fprintf('  PASS: RayTracePropagator dz stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator dz stored\n');
    failed = failed + 1;
end

% testRayTraceZ0Stored
rp_z0 = RayTracePropagator(grid, 'RK4', [], 0.01);
if (rp_z0.z0 == 0.01)
    fprintf('  PASS: RayTracePropagator z0 stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayTracePropagator z0 stored\n');
    failed = failed + 1;
end

% -----------------------------------------------------------------
% Polymorphism test — swap propagators without changing caller
% -----------------------------------------------------------------

% testPolymorphicSwap
propagators = {AnalyticPropagator(grid), FFTPropagator(grid, lambda)};
all_ok = true;
for i = 1:numel(propagators)
    f = propagators{i}.propagate(beam, 0.05);
    if ~(size(f, 1) == 64 && all(all(isfinite(f))))
        all_ok = false;
    end
end
if all_ok
    fprintf('  PASS: polymorphic propagate (analytic + FFT)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: polymorphic propagate\n');
    failed = failed + 1;
end

fprintf('\n=== Propagators: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
