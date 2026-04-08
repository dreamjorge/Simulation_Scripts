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

% testHermitePolyN0
x_test = linspace(-2, 2, 10);
H0 = HermiteBeam.hermitePoly(0, x_test);
if (all(abs(H0 - 1) < 1e-10))
    fprintf('  PASS: hermitePoly n=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly n=0\n');
    failed = failed + 1;
end

% testHermitePolyN1
H1 = HermiteBeam.hermitePoly(1, x_test);
expected_H1 = 2 * x_test;
if (all(abs(H1 - expected_H1) < 1e-10))
    fprintf('  PASS: hermitePoly n=1\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly n=1\n');
    failed = failed + 1;
end

% testHermitePolyN2
H2 = HermiteBeam.hermitePoly(2, x_test);
expected_H2 = 4 * x_test.^2 - 2;
if (all(abs(H2 - expected_H2) < 1e-10))
    fprintf('  PASS: hermitePoly n=2\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly n=2\n');
    failed = failed + 1;
end

% testHermitePolyN3
H3 = HermiteBeam.hermitePoly(3, x_test);
expected_H3 = 8 * x_test.^3 - 12 * x_test;
if (all(abs(H3 - expected_H3) < 1e-10))
    fprintf('  PASS: hermitePoly n=3\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly n=3\n');
    failed = failed + 1;
end

% testHermitePolyMatrixInput
x_mat = rand(5, 5);
H_mat = HermiteBeam.hermitePoly(2, x_mat);
if (size(H_mat) == size(x_mat))
    fprintf('  PASS: hermitePoly matrix input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly matrix input\n');
    failed = failed + 1;
end

% testPhaseIncluded
hp_phase = HermiteParameters(0.05, w0, lambda, 1, 1);
hb_phase = HermiteBeam(X, Y, hp_phase);
field_complex = hb_phase.OpticalField;
if (~isreal(field_complex) || any(any(abs(imag(field_complex)) > 0)))
    fprintf('  PASS: phase included in field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase included in field\n');
    failed = failed + 1;
end

% testHermitePolyN4
H4 = HermiteBeam.hermitePoly(4, x_test);
expected_H4 = 16 * x_test.^4 - 48 * x_test.^2 + 12;
if (all(abs(H4 - expected_H4) < 1e-10))
    fprintf('  PASS: hermitePoly n=4\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly n=4\n');
    failed = failed + 1;
end

% testFieldAtDifferentZ
hp_z1 = HermiteParameters(0.01, w0, lambda, 1, 1);
hp_z2 = HermiteParameters(0.1, w0, lambda, 1, 1);
hb_z1 = HermiteBeam(X, Y, hp_z1);
hb_z2 = HermiteBeam(X, Y, hp_z2);
if (all(all(isfinite(hb_z1.OpticalField))) && all(all(isfinite(hb_z2.OpticalField))))
    fprintf('  PASS: field at different z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at different z\n');
    failed = failed + 1;
end

% testHigherNMOrders
hb_hn = HermiteBeam(X, Y, HermiteParameters(0, w0, lambda, 3, 3));
if (all(all(isfinite(hb_hn.OpticalField))))
    fprintf('  PASS: higher n m orders\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher n m\n');
    failed = failed + 1;
end

% testZeroOrderPhaseOnly
hb_zp = HermiteBeam(X, Y, HermiteParameters(0, w0, lambda, 0, 0));
if (hb_zp.Parameters.PhiPhase == 0)
    fprintf('  PASS: zero order phase only\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order phase\n');
    failed = failed + 1;
end

% testCoordinatesNotEmpty
if (~isempty(hb.x) && ~isempty(hb.y))
    fprintf('  PASS: coordinates not empty\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates empty\n');
    failed = failed + 1;
end

% testFieldAmplitude
hb_amp = HermiteBeam(X, Y, HermiteParameters(0, w0, lambda, 0, 0));
max_amp = max(max(abs(hb_amp.OpticalField)));
if (max_amp > 0)
    fprintf('  PASS: field amplitude positive\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field amplitude\n');
    failed = failed + 1;
end

% testDifferentWavelengths
hb_l1 = HermiteBeam(X, Y, HermiteParameters(0, w0, 532e-9, 1, 1));
hb_l2 = HermiteBeam(X, Y, HermiteParameters(0, w0, 1064e-9, 1, 1));
if (all(all(isfinite(hb_l1.OpticalField))) && all(all(isfinite(hb_l2.OpticalField))))
    fprintf('  PASS: different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different wavelengths\n');
    failed = failed + 1;
end

% testHermitePolyNegativeValues
x_neg = linspace(-1, 1, 5);
H_neg = HermiteBeam.hermitePoly(2, x_neg);
if (all(isfinite(H_neg)))
    fprintf('  PASS: hermitePoly negative values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: hermitePoly negative\n');
    failed = failed + 1;
end

% testAsymmetricOrders
hb_asym = HermiteBeam(X, Y, HermiteParameters(0, w0, lambda, 2, 1));
if (all(all(isfinite(hb_asym.OpticalField))))
    fprintf('  PASS: asymmetric orders\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric orders\n');
    failed = failed + 1;
end

% testWaistFromParameters
if (hb.Parameters.Waist > 0)
    fprintf('  PASS: waist from parameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist from params\n');
    failed = failed + 1;
end

% testGouyPhaseIncluded
hp_gouy = HermiteParameters(0.05, w0, lambda, 1, 1);
hb_gouy = HermiteBeam(X, Y, hp_gouy);
if (hb_gouy.Parameters.GouyPhase > 0)
    fprintf('  PASS: Gouy phase included\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Gouy phase\n');
    failed = failed + 1;
end

fprintf('\n=== HermiteBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
