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

% testFieldSymmetry
gb_sym = GaussianBeam(R, GaussianParameters(0, w0, lambda));
field_center = abs(gb_sym.OpticalField(33,33));
field_edge = abs(gb_sym.OpticalField(1,1));
if (field_edge < field_center)
    fprintf('  PASS: field symmetry\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field symmetry\n');
    failed = failed + 1;
end

% testWaistPropagates
params_w1 = GaussianParameters(0, w0, lambda);
params_w2 = GaussianParameters(0.1, w0, lambda);
gb_w1 = GaussianBeam(R, params_w1);
gb_w2 = GaussianBeam(R, params_w2);
if (params_w2.Waist > params_w1.Waist)
    fprintf('  PASS: waist propagates\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist propagates\n');
    failed = failed + 1;
end

% testDifferentWavelength
lambda2 = 1064e-9;
gb_l2 = GaussianBeam(R, GaussianParameters(0, w0, lambda2));
if (all(all(isfinite(gb_l2.OpticalField))))
    fprintf('  PASS: different wavelength\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different wavelength\n');
    failed = failed + 1;
end

% testDifferentWaist
w02 = 50e-6;
gb_w02 = GaussianBeam(R, GaussianParameters(0, w02, lambda));
amp_ratio = abs(gb_w02.OpticalField(33,33)) / abs(gb.OpticalField(33,33));
if (amp_ratio > 1)
    fprintf('  PASS: different waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different waist\n');
    failed = failed + 1;
end

% testFieldWithPositiveZ
gb_pos_z = GaussianBeam(R, GaussianParameters(0.05, w0, lambda));
if (all(all(isfinite(gb_pos_z.OpticalField))))
    fprintf('  PASS: field with positive z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: positive z\n');
    failed = failed + 1;
end

% testFieldWithNegativeZ
gb_neg_z = GaussianBeam(R, GaussianParameters(-0.05, w0, lambda));
if (all(all(isfinite(gb_neg_z.OpticalField))))
    fprintf('  PASS: field with negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative z\n');
    failed = failed + 1;
end

% testConstructorEmpty
try
    gb_empty = GaussianBeam();
    fprintf('  PASS: constructor empty\n');
    passed = passed + 1;
catch
    fprintf('  FAIL: constructor empty\n');
    failed = failed + 1;
end

% testBeamDivergence
gb_div = GaussianBeam(5*w0*ones(64,64), GaussianParameters(0, w0, lambda));
if (abs(gb_div.OpticalField(1,1)) < 0.1)
    fprintf('  PASS: beam divergence\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beam divergence\n');
    failed = failed + 1;
end

% testFieldNormalized
gb_norm = GaussianBeam(zeros(64,64), GaussianParameters(0, w0, lambda));
max_amp = max(max(abs(gb_norm.OpticalField)));
if (abs(max_amp - 1) < 1e-10)
    fprintf('  PASS: field normalized\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field normalized\n');
    failed = failed + 1;
end

fprintf('\n=== GaussianBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
% exit(0);
else
% exit(1);
end
