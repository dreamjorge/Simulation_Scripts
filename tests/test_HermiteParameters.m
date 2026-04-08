#!/usr/bin/env octave
% Tests for HermiteParameters

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== HermiteParameters Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
z = 0.05;

% testConstructor
hp = HermiteParameters(z, w0, lambda, 1, 1);
if (hp.n == 1 && hp.m == 1)
    fprintf('  PASS: constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor\n');
    failed = failed + 1;
end

% testDefaultNandM
hp0 = HermiteParameters(z, w0, lambda);
if (hp0.n == 0 && hp0.m == 0)
    fprintf('  PASS: default n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: default n m\n');
    failed = failed + 1;
end

% testHermiteWaistX
hp1 = HermiteParameters(z, w0, lambda, 2, 0);
expected_wx = hp1.Waist * sqrt(2 + 1);
if (abs(hp1.HermiteWaistX - expected_wx) / expected_wx < 1e-10)
    fprintf('  PASS: HermiteWaistX\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HermiteWaistX\n');
    failed = failed + 1;
end

% testHermiteWaistY
expected_wy = hp1.Waist * sqrt(0 + 1);
if (abs(hp1.HermiteWaistY - expected_wy) / expected_wy < 1e-10)
    fprintf('  PASS: HermiteWaistY\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HermiteWaistY\n');
    failed = failed + 1;
end

% testHermiteWaist
hp2 = HermiteParameters(z, w0, lambda, 2, 3);
expected_w = hp2.Waist * sqrt(2 + 3 + 1);
if (abs(hp2.HermiteWaist - expected_w) / expected_w < 1e-10)
    fprintf('  PASS: HermiteWaist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HermiteWaist\n');
    failed = failed + 1;
end

% testPhiPhase
zr = hp2.RayleighDistance;
expected_phi = (2 + 3) * atan(z/zr);
if (abs(hp2.PhiPhase - expected_phi) < 1e-10)
    fprintf('  PASS: PhiPhase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase\n');
    failed = failed + 1;
end

% testPhiPhaseZeroAtWaist
hp_z0 = HermiteParameters(0, w0, lambda, 1, 1);
if (abs(hp_z0.PhiPhase) < 1e-10)
    fprintf('  PASS: PhiPhase zero at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase at waist\n');
    failed = failed + 1;
end

% testStaticGetWaistOneDirection
zr_s = pi*w0^2/lambda;
wH_s = HermiteParameters.getWaistOneDirection(0.05, w0, zr_s, 1);
if (wH_s > 0)
    fprintf('  PASS: static getWaistOneDirection\n');
    passed = passed + 1;
else
    fprintf('  FAIL: static getWaistOneDirection\n');
    failed = failed + 1;
end

% testStaticGetWaist
wH2_s = HermiteParameters.getWaist(0.05, w0, zr_s, 1, 2);
if (wH2_s > 0)
    fprintf('  PASS: static getWaist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: static getWaist\n');
    failed = failed + 1;
end

% testWaistIncreasesWithZ
hp_z1 = HermiteParameters(0.01, w0, lambda, 1, 1);
hp_z2 = HermiteParameters(0.1, w0, lambda, 1, 1);
if (hp_z2.HermiteWaist > hp_z1.HermiteWaist)
    fprintf('  PASS: waist grows with z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist grows with z\n');
    failed = failed + 1;
end

% testHigherOrdersLargerWaist
hp_12 = HermiteParameters(z, w0, lambda, 1, 2);
hp_34 = HermiteParameters(z, w0, lambda, 3, 4);
if (hp_34.HermiteWaist > hp_12.HermiteWaist)
    fprintf('  PASS: higher orders larger waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher orders\n');
    failed = failed + 1;
end

% testWavelengthProperty
if (hp.Wavelength == lambda)
    fprintf('  PASS: Wavelength stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Wavelength\n');
    failed = failed + 1;
end

fprintf('\n=== HermiteParameters: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
