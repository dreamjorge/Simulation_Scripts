#!/usr/bin/env octave
% Tests for GaussianBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== GaussianBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
[R, ~] = cart2pol(X, Y);

% testFieldGeneration
params = GaussianParameters(0.01, w0, lambda);
gb = GaussianBeam(R, params);
if (size(gb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAmplitudeAtCenter
params0 = GaussianParameters(0, w0, lambda);
gb0 = GaussianBeam(zeros(64,64), params0);
expected_amp = 1;
if (abs(abs(gb0.OpticalField(33,33)) - expected_amp) < 1e-10)
    fprintf('  PASS: amplitude at center\n');
    passed = passed + 1;
else
    fprintf('  FAIL: amplitude at center\n');
    failed = failed + 1;
end

% testDecaysWithRadius
R_far = 3*w0*ones(64,64);
gb_far = GaussianBeam(R_far, params0);
if (abs(gb_far.OpticalField(1,1)) < 0.01)
    fprintf('  PASS: decays with radius\n');
    passed = passed + 1;
else
    fprintf('  FAIL: decays with radius\n');
    failed = failed + 1;
end

% testAtZ0Finite
params_z0 = GaussianParameters(0, w0, lambda);
gb_z0 = GaussianBeam(R, params_z0);
if (all(all(isfinite(gb_z0.OpticalField))))
    fprintf('  PASS: at z=0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z=0 finite\n');
    failed = failed + 1;
end

% testParametersStored
if (gb.Parameters.InitialWaist == w0)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters stored\n');
    failed = failed + 1;
end

% testComplexField
params_c = GaussianParameters(0.05, w0, lambda);
gb_c = GaussianBeam(R, params_c);
if (~isreal(gb_c.OpticalField))
    fprintf('  PASS: complex field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: complex field\n');
    failed = failed + 1;
end

% testPhaseAtWaist
gb_p = GaussianBeam(zeros(64,64), GaussianParameters(0, w0, lambda));
phase_center = angle(gb_p.OpticalField(33,33));
if (abs(phase_center) < pi/2)
    fprintf('  PASS: phase at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase at waist\n');
    failed = failed + 1;
end

% testzCoordinateStored
gb_z = GaussianBeam(R, GaussianParameters(0.1, w0, lambda));
if (gb_z.Parameters.zCoordinate == 0.1)
    fprintf('  PASS: zCoordinate stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zCoordinate stored\n');
    failed = failed + 1;
end

fprintf('\n=== GaussianBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
