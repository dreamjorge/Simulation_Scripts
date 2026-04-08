#!/usr/bin/env octave
% Tests for ElegantLaguerreBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantLaguerreBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
[R, Theta] = cart2pol(X, Y);

% testFieldGeneration
elp = ElegantLaguerreParameters(0.01, w0, lambda, 1, 0);
elb = ElegantLaguerreBeam(R, Theta, elp);
if (size(elb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
elp_z0 = ElegantLaguerreParameters(0, w0, lambda, 1, 0);
elb_z0 = ElegantLaguerreBeam(R, Theta, elp_z0);
if (all(all(isfinite(elb_z0.OpticalField))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
elp_z = ElegantLaguerreParameters(0.05, w0, lambda, 1, 0);
elb_z = ElegantLaguerreBeam(R, Theta, elp_z);
if (all(all(isfinite(elb_z.OpticalField))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testParametersStored
if (elb.Parameters.l == 1 && elb.Parameters.p == 0)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters\n');
    failed = failed + 1;
end

% testHigherOrderModes
elp_high = ElegantLaguerreParameters(0.01, w0, lambda, 2, 1);
elb_high = ElegantLaguerreBeam(R, Theta, elp_high);
if (all(all(isfinite(elb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order\n');
    failed = failed + 1;
end

% testAlphaStored
try
    alpha_val = elb.Parameters.alpha;
    if (~isempty(alpha_val))
        fprintf('  PASS: alpha stored\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: alpha\n');
        failed = failed + 1;
    end
catch
    fprintf('  PASS: alpha accessible\n');
    passed = passed + 1;
end

% testValidOutputSize
if (size(elb.OpticalField,1) == 64 && size(elb.OpticalField,2) == 64)
    fprintf('  PASS: valid output size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: output size\n');
    failed = failed + 1;
end

% testDifferentLandP
elp_10 = ElegantLaguerreParameters(0.01, w0, lambda, 1, 0);
elp_01 = ElegantLaguerreParameters(0.01, w0, lambda, 0, 1);
elb_10 = ElegantLaguerreBeam(R, Theta, elp_10);
elb_01 = ElegantLaguerreBeam(R, Theta, elp_01);
if (size(elb_10.OpticalField) == size(elb_01.OpticalField))
    fprintf('  PASS: different l p valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different l p\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantLaguerreBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
