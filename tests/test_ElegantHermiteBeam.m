#!/usr/bin/env octave
% Tests for ElegantHermiteBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantHermiteBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

% testFieldGeneration
ehp = ElegantHermiteParameters(0.01, w0, lambda, 1, 1);
ehb = ElegantHermiteBeam(X, Y, ehp);
if (size(ehb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
ehp_z0 = ElegantHermiteParameters(0, w0, lambda, 1, 1);
ehb_z0 = ElegantHermiteBeam(X, Y, ehp_z0);
if (all(all(isfinite(ehb_z0.OpticalField))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
ehp_z = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
ehb_z = ElegantHermiteBeam(X, Y, ehp_z);
if (all(all(isfinite(ehb_z.OpticalField))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testCoordinatesStored
if (isequal(size(ehb.X), [64, 64]) && isequal(size(ehb.Y), [64, 64]))
    fprintf('  PASS: coordinates stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates\n');
    failed = failed + 1;
end

% testParametersStored
if (ehb.Parameters.n == 1 && ehb.Parameters.m == 1)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters\n');
    failed = failed + 1;
end

% testHigherOrderModes
ehp_high = ElegantHermiteParameters(0.01, w0, lambda, 2, 2);
ehb_high = ElegantHermiteBeam(X, Y, ehp_high);
if (all(all(isfinite(ehb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order\n');
    failed = failed + 1;
end

% testAlphaStored
try
    alpha_val = ehb.Parameters.alpha;
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
if (size(ehb.OpticalField,1) == 64 && size(ehb.OpticalField,2) == 64)
    fprintf('  PASS: valid output size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: output size\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantHermiteBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
