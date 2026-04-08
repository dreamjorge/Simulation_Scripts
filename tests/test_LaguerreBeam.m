#!/usr/bin/env octave
% Tests for LaguerreBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== LaguerreBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
[R, Theta] = cart2pol(X, Y);

% testFieldGeneration
lp = LaguerreParameters(0, w0, lambda, 1, 0);
lb = LaguerreBeam(R, Theta, lp);
if (size(lb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testHigherOrderModes
lp_high = LaguerreParameters(0, w0, lambda, 2, 1);
lb_high = LaguerreBeam(R, Theta, lp_high);
if (all(all(isfinite(lb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order modes\n');
    failed = failed + 1;
end

% testCoordinatesStored
if (isequal(size(lb.r), [64, 64]) && isequal(size(lb.theta), [64, 64]))
    fprintf('  PASS: coordinates stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates stored\n');
    failed = failed + 1;
end

% testParametersStored
if (lb.Parameters.l == 1 && lb.Parameters.p == 0)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters stored\n');
    failed = failed + 1;
end

% testValidFieldAtWaist
lp_z0 = LaguerreParameters(0, w0, lambda, 1, 0);
lb_z0 = LaguerreBeam(R, Theta, lp_z0);
if (all(all(isfinite(lb_z0.OpticalField))))
    fprintf('  PASS: valid field at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: valid field at waist\n');
    failed = failed + 1;
end

% testDifferentLandP
lp_10 = LaguerreParameters(0, w0, lambda, 1, 0);
lp_01 = LaguerreParameters(0, w0, lambda, 0, 1);
lb_10 = LaguerreBeam(R, Theta, lp_10);
lb_01 = LaguerreBeam(R, Theta, lp_01);
if (size(lb_10.OpticalField) == size(lb_01.OpticalField))
    fprintf('  PASS: different l p valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different l p\n');
    failed = failed + 1;
end

% testLaguerreZeroPZero
lp_00 = LaguerreParameters(0, w0, lambda, 0, 0);
lb_00 = LaguerreBeam(R, Theta, lp_00);
gb = GaussianBeam(R, GaussianParameters(0, w0, lambda));
if (size(lb_00.OpticalField) == size(gb.OpticalField))
    fprintf('  PASS: l=0 p=0 valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: l=0 p=0\n');
    failed = failed + 1;
end

fprintf('\n=== LaguerreBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
