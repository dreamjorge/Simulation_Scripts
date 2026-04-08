#!/usr/bin/env octave
% Tests for LaguerreParameters

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== LaguerreParameters Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
z = 0.05;

% testConstructor
lp = LaguerreParameters(z, w0, lambda, 1, 1);
if (lp.l == 1 && lp.p == 1)
    fprintf('  PASS: constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor\n');
    failed = failed + 1;
end

% testDefaultLandP
lp0 = LaguerreParameters(z, w0, lambda);
if (lp0.l == 0 && lp0.p == 0)
    fprintf('  PASS: default l p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: default l p\n');
    failed = failed + 1;
end

% testLaguerreWaist
lp1 = LaguerreParameters(z, w0, lambda, 1, 1);
expected_lw = lp1.Waist * sqrt(2*1 + abs(1) + 1);
if (abs(lp1.LaguerreWaist - expected_lw) / expected_lw < 1e-10)
    fprintf('  PASS: LaguerreWaist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: LaguerreWaist\n');
    failed = failed + 1;
end

% testPhiPhase
zr = lp1.RayleighDistance;
expected_phi = (abs(1) + 2*1) * atan(z/zr);
if (abs(lp1.PhiPhase - expected_phi) < 1e-10)
    fprintf('  PASS: PhiPhase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase\n');
    failed = failed + 1;
end

% testPhiPhaseZeroAtWaist
lp_z0 = LaguerreParameters(0, w0, lambda, 1, 1);
if (abs(lp_z0.PhiPhase) < 1e-10)
    fprintf('  PASS: PhiPhase zero at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase at waist\n');
    failed = failed + 1;
end

% testStaticGetWaist
zr_s = pi*w0^2/lambda;
wL_s = LaguerreParameters.getWaist(0.05, w0, zr_s, 1, 1);
if (wL_s > 0)
    fprintf('  PASS: static getWaist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: static getWaist\n');
    failed = failed + 1;
end

% testAssociatedLaguerrePolynomial
x_test = [0, 0.5, 1, 2];
L_test = LaguerreParameters.getAssociatedLaguerrePolynomial(2, 1, x_test);
if (numel(L_test) == 4 && all(isfinite(L_test)))
    fprintf('  PASS: getAssociatedLaguerrePolynomial\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getAssociatedLaguerrePolynomial\n');
    failed = failed + 1;
end

% testNegativeL
lp_neg = LaguerreParameters(z, w0, lambda, -1, 1);
if (lp_neg.l == -1 && lp_neg.p == 1)
    fprintf('  PASS: negative l\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l\n');
    failed = failed + 1;
end

% testWaistIncreasesWithZ
lp_z1 = LaguerreParameters(0.01, w0, lambda, 1, 0);
lp_z2 = LaguerreParameters(0.1, w0, lambda, 1, 0);
if (lp_z2.LaguerreWaist > lp_z1.LaguerreWaist)
    fprintf('  PASS: waist grows with z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist grows with z\n');
    failed = failed + 1;
end

% testHigherLLargerWaist
lp_l1 = LaguerreParameters(z, w0, lambda, 1, 0);
lp_l2 = LaguerreParameters(z, w0, lambda, 2, 0);
if (lp_l2.LaguerreWaist > lp_l1.LaguerreWaist)
    fprintf('  PASS: higher l larger waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher l\n');
    failed = failed + 1;
end

% testWavelengthProperty
if (lp.Wavelength == lambda)
    fprintf('  PASS: Wavelength stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Wavelength\n');
    failed = failed + 1;
end

fprintf('\n=== LaguerreParameters: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
