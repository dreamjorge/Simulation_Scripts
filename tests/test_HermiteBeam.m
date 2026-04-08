#!/usr/bin/env octave
% Tests for HermiteBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== HermiteBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

% testFieldGeneration
hp = HermiteParameters(0, w0, lambda, 1, 1);
hb = HermiteBeam(X, Y, hp);
if (size(hb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testHigherOrderModes
hp_high = HermiteParameters(0, w0, lambda, 2, 3);
hb_high = HermiteBeam(X, Y, hp_high);
if (all(all(isfinite(hb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order modes\n');
    failed = failed + 1;
end

% testZeroOrderEqualsGaussian
[R, ~] = cart2pol(X, Y);
hp_zero = HermiteParameters(0, w0, lambda, 0, 0);
hb_zero = HermiteBeam(X, Y, hp_zero);
gb = GaussianBeam(R, GaussianParameters(0, w0, lambda));
diff_center = abs(hb_zero.OpticalField(33,33) - gb.OpticalField(33,33));
if (diff_center < 1e-10)
    fprintf('  PASS: zero order equals Gaussian\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order equals Gaussian\n');
    failed = failed + 1;
end

% testCoordinatesStored
if (isequal(size(hb.x), [64, 64]) && isequal(size(hb.y), [64, 64]))
    fprintf('  PASS: coordinates stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates stored\n');
    failed = failed + 1;
end

% testParametersStored
if (hb.Parameters.n == 1 && hb.Parameters.m == 1)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters stored\n');
    failed = failed + 1;
end

% testValidFieldAtWaist
hp_z0 = HermiteParameters(0, w0, lambda, 1, 1);
hb_z0 = HermiteBeam(X, Y, hp_z0);
if (all(all(isfinite(hb_z0.OpticalField))))
    fprintf('  PASS: valid field at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: valid field at waist\n');
    failed = failed + 1;
end

% testDifferentOrders
hp_10 = HermiteParameters(0, w0, lambda, 1, 0);
hp_01 = HermiteParameters(0, w0, lambda, 0, 1);
hb_10 = HermiteBeam(X, Y, hp_10);
hb_01 = HermiteBeam(X, Y, hp_01);
if (size(hb_10.OpticalField) == size(hb_01.OpticalField))
    fprintf('  PASS: different orders valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different orders\n');
    failed = failed + 1;
end

fprintf('\n=== HermiteBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
