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

% testNegativeL
lp_neg = LaguerreParameters(0, w0, lambda, -1, 0);
lb_neg = LaguerreBeam(R, Theta, lp_neg);
if (all(all(isfinite(lb_neg.OpticalField))))
    fprintf('  PASS: negative l valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l\n');
    failed = failed + 1;
end

% testHigherP
lp_p2 = LaguerreParameters(0, w0, lambda, 0, 2);
lb_p2 = LaguerreBeam(R, Theta, lp_p2);
if (all(all(isfinite(lb_p2.OpticalField))))
    fprintf('  PASS: higher p order valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher p order\n');
    failed = failed + 1;
end

% testPhaseVariation
lp_phase = LaguerreParameters(0.05, w0, lambda, 1, 0);
lb_phase = LaguerreBeam(R, Theta, lp_phase);
if (lp_phase.PhiPhase ~= 0)
    fprintf('  PASS: phase variation included\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase variation\n');
    failed = failed + 1;
end

% testAzimuthalPhase
lp_az = LaguerreParameters(0, w0, lambda, 1, 0);
lb_az = LaguerreBeam(R, Theta, lp_az);
theta_sample = Theta(32, 32);
phase_at_theta = angle(lb_az.OpticalField(32,32));
expected_phase_contrib = 1 * theta_sample;
if (abs(phase_at_theta - expected_phase_contrib) < pi)
    fprintf('  PASS: azimuthal phase contribution\n');
    passed = passed + 1;
else
    fprintf('  FAIL: azimuthal phase\n');
    failed = failed + 1;
end

% testAtPropagation
lp_prop = LaguerreParameters(0.1, w0, lambda, 1, 0);
lb_prop = LaguerreBeam(R, Theta, lp_prop);
if (all(all(isfinite(lb_prop.OpticalField))))
    fprintf('  PASS: field at propagation valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at propagation\n');
    failed = failed + 1;
end

% testCombinedLandP
lp_comb = LaguerreParameters(0, w0, lambda, 2, 3);
lb_comb = LaguerreBeam(R, Theta, lp_comb);
if (all(all(isfinite(lb_comb.OpticalField))))
    fprintf('  PASS: combined l p valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: combined l p\n');
    failed = failed + 1;
end

% testFieldAmplitude
lb_amp = LaguerreBeam(R, Theta, LaguerreParameters(0, w0, lambda, 0, 0));
max_amp = max(max(abs(lb_amp.OpticalField)));
if (max_amp > 0)
    fprintf('  PASS: field amplitude positive\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field amplitude\n');
    failed = failed + 1;
end

% testWaistFromParameters
if (lb.Parameters.Waist > 0)
    fprintf('  PASS: waist from parameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist from params\n');
    failed = failed + 1;
end

% testLaguerreWaist
lp_lw = LaguerreParameters(0, w0, lambda, 1, 1);
lb_lw = LaguerreBeam(R, Theta, lp_lw);
if (lb_lw.Parameters.LaguerreWaist > 0)
    fprintf('  PASS: Laguerre waist valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Laguerre waist\n');
    failed = failed + 1;
end

% testDifferentWavelengths
lb_l1 = LaguerreBeam(R, Theta, LaguerreParameters(0, w0, 532e-9, 1, 0));
lb_l2 = LaguerreBeam(R, Theta, LaguerreParameters(0, w0, 1064e-9, 1, 0));
if (all(all(isfinite(lb_l1.OpticalField))) && all(all(isfinite(lb_l2.OpticalField))))
    fprintf('  PASS: different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different wavelengths\n');
    failed = failed + 1;
end

% testCoordinatesNotEmpty
if (~isempty(lb.r) && ~isempty(lb.theta))
    fprintf('  PASS: coordinates not empty\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates empty\n');
    failed = failed + 1;
end

% testFieldAtDifferentZ
lp_z1 = LaguerreParameters(0.01, w0, lambda, 1, 0);
lp_z2 = LaguerreParameters(0.1, w0, lambda, 1, 0);
lb_z1 = LaguerreBeam(R, Theta, lp_z1);
lb_z2 = LaguerreBeam(R, Theta, lp_z2);
if (all(all(isfinite(lb_z1.OpticalField))) && all(all(isfinite(lb_z2.OpticalField))))
    fprintf('  PASS: field at different z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at different z\n');
    failed = failed + 1;
end

% testNegativeLWithP
lp_negp = LaguerreParameters(0, w0, lambda, -2, 1);
lb_negp = LaguerreBeam(R, Theta, lp_negp);
if (all(all(isfinite(lb_negp.OpticalField))))
    fprintf('  PASS: negative l with p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l with p\n');
    failed = failed + 1;
end

% testPhiPhaseAtWaist
lp_pw = LaguerreParameters(0, w0, lambda, 1, 1);
if (abs(lp_pw.PhiPhase) < 1e-10)
    fprintf('  PASS: PhiPhase at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase at waist\n');
    failed = failed + 1;
end

% testHigherLLargerWaist
lp_ll = LaguerreParameters(0, w0, lambda, 2, 0);
lb_ll = LaguerreBeam(R, Theta, lp_ll);
if (lb_ll.Parameters.LaguerreWaist > lb.Parameters.LaguerreWaist)
    fprintf('  PASS: higher l larger waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher l waist\n');
    failed = failed + 1;
end

fprintf('\n=== LaguerreBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
