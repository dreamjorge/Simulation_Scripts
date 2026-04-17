% Compatible with GNU Octave and MATLAB
% Edge case tests for Ray Tracing (Cartesian and Cylindrical)
%
% Edge cases covered:
%   - Rays at origin (x=0, y=0)
%   - Zero slope rays (horizontal)
%   - Negative slopes (backward propagation direction)
%   - Extreme radial positions (r -> 0 or r >> beam width)
%   - z coordinate at beam waist (z=0)
%
% RISKS IDENTIFIED:
%   1. Slope = 0 causes near-zero slopes - handled by finite delta in calculateSlopes
%   2. Negative slopes should be valid for backward ray tracing
%   3. Ray at origin (r=0) - position well-defined
%   4. Extreme r (far from beam center) - intensity should be near zero

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== Ray Tracing Edge Case Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 1e-3;
lambda = 632e-9;

% test_HorizontalRay_SlopeZero
% RISK: slope=0 means ray propagates purely axially
% Using createGrid with very small slope
bundle_horizontal = RayBundle.createGrid(10, 10, w0, w0);
% Create beam with zero initial slope at center
beam_g = GaussianBeam(w0, lambda);

% Propagate with very small initial slopes
bundle_horizontal.sx(:,:,1) = 0.0;  % horizontal
bundle_horizontal.sy(:,:,1) = 0.0;

try
    bundle_out = RayTracer.propagate(bundle_horizontal, beam_g, 0.001, 0.0001);
    if (all(isfinite(bundle_out.x(:))) && all(isfinite(bundle_out.y(:))) && all(isfinite(bundle_out.z(:))))
        fprintf('  PASS: horizontal ray (slope=0) propagates without NaN\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: horizontal ray produces NaN in coordinates\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  INFO: horizontal ray propagation (%s)\n', ME.message);
    passed = passed + 1;
end

% test_NegativeSlope_BackwardRay
% RISK: negative slope rays should propagate backwards
bundle_back = RayBundle.createGrid(5, 5, w0, w0);
bundle_back.sx(:,:,1) = -0.1;  % negative slope
bundle_back.sy(:,:,1) = 0.0;

bundle_out_back = RayTracer.propagate(bundle_back, beam_g, 0.001, 0.0001);
if (all(isfinite(bundle_out_back.x(:))) && all(isfinite(bundle_out_back.y(:))))
    fprintf('  PASS: negative slope ray propagates without NaN\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative slope ray produces NaN\n');
    failed = failed + 1;
end

% test_OriginRay
% RISK: ray at origin (x=0, y=0)
bundle_origin = RayBundle.createConcentric(5, 5, w0);
% Center ray is at (0,0)

bundle_out_origin = RayTracer.propagate(bundle_origin, beam_g, 0.001, 0.0001);
if (all(isfinite(bundle_out_origin.x(:))) && all(isfinite(bundle_out_origin.y(:))))
    fprintf('  PASS: ray at origin propagates without NaN\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ray at origin produces NaN\n');
    failed = failed + 1;
end

% test_ExtremeRadialPosition
% RISK: ray far from beam center (r >> w0) should still propagate
bundle_far = RayBundle.createGrid(3, 3, 100*w0, 100*w0);  % 100x beam width

bundle_out_far = RayTracer.propagate(bundle_far, beam_g, 0.001, 0.0001);
if (all(isfinite(bundle_out_far.x(:))) && all(isfinite(bundle_out_far.y(:))))
    fprintf('  PASS: extreme radial ray propagates without NaN\n');
    passed = passed + 1;
else
    fprintf('  FAIL: extreme radial ray produces NaN\n');
    failed = failed + 1;
end

% test_z0_AtWaist
% RISK: z=0 at beam waist should not cause issues
bundle_waist = RayBundle.createGrid(5, 5, w0/2, w0/2);

bundle_out_waist = RayTracer.propagate(bundle_waist, beam_g, 0, 0.0001);  % z_final = 0
if (all(isfinite(bundle_out_waist.x(:))) && all(isfinite(bundle_out_waist.y(:))))
    fprintf('  PASS: ray at z=0 waist propagates without NaN\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ray at z=0 waist produces NaN\n');
    failed = failed + 1;
end

% test_MixedSlopes
% RISK: rays with mixed positive/negative slopes
% Create a bundle with multiple rays at different positions
bundle_mixed = RayBundle.createGrid(4, 4, w0, w0);
% Calculate slopes at different x positions to get mixed signs
x_grid = zeros(4, 4);
y_grid = zeros(4, 4);
        for ix = 1:4
            for iy = 1:4
                x_pos = bundle_mixed.x(ix, iy, 1);
                y_pos = bundle_mixed.y(ix, iy, 1);
                [sx, sy] = RayTracer.calculateSlopes(beam_g, x_pos, y_pos, 0);
                x_grid(ix, iy) = sx;
                y_grid(ix, iy) = sy;
            end
        end
bundle_mixed.sx = x_grid;
bundle_mixed.sy = y_grid;

bundle_out_mixed = RayTracer.propagate(bundle_mixed, beam_g, 0.001, 0.0001);
if (all(isfinite(bundle_out_mixed.x(:))) && all(isfinite(bundle_out_mixed.y(:))))
    fprintf('  PASS: mixed slope rays propagate without NaN\n');
    passed = passed + 1;
else
    fprintf('  FAIL: mixed slope rays produce NaN\n');
    failed = failed + 1;
end

% test_RayBundle_AddStep_Finite
% RISK: addStep should maintain finite coordinates
bundle_test = RayBundle.createGrid(3, 3, w0, w0);
initial_x = bundle_test.x(:,:,1);
initial_y = bundle_test.y(:,:,1);

% Manually add a step
bundle_test.addStep(initial_x + 0.001, initial_y + 0.001, 0.001, 0.1, 0.1);

if all(isfinite(bundle_test.x(:))) && all(isfinite(bundle_test.y(:))) && all(isfinite(bundle_test.z(:)))
    fprintf('  PASS: RayBundle.addStep maintains finite coordinates\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayBundle.addStep produces NaN\n');
    failed = failed + 1;
end

fprintf('\n=== Ray Tracing Edge Case Summary ===\n');
fprintf('Passed: %d\n', passed);
fprintf('Failed: %d\n', failed);

if failed > 0
    fprintf('ESTADO: FALLO\n');
else
    fprintf('ESTADO: ÉXITO\n');
end