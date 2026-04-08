#!/usr/bin/env octave
% Tests for PhysicalConstants

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== PhysicalConstants Tests ===\n\n');
passed = 0;
failed = 0;

lambda = 632.8e-9;
w0 = 100e-6;
zr = pi * w0^2 / lambda;

% testWaveNumber
k = PhysicalConstants.waveNumber(lambda);
expected = 2*pi / lambda;
if (abs(k - expected) < 1e-5)
    fprintf('  PASS: waveNumber\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber\n');
    failed = failed + 1;
end

% testRayleighDistance
zr_test = PhysicalConstants.rayleighDistance(w0, lambda);
expected_zr = pi * w0^2 / lambda;
if (abs(zr_test - expected_zr) < 1e-10)
    fprintf('  PASS: rayleighDistance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: rayleighDistance\n');
    failed = failed + 1;
end

% testWaistAtZ
w = PhysicalConstants.waistAtZ(w0, 0.05, lambda);
expected_w = w0 * sqrt(1 + (0.05/zr)^2);
if (abs(w - expected_w) / expected_w < 1e-5)
    fprintf('  PASS: waistAtZ\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ\n');
    failed = failed + 1;
end

% testRadiusOfCurvature
R = PhysicalConstants.radiusOfCurvature(0.05, zr);
if (R > 0)
    fprintf('  PASS: radiusOfCurvature\n');
    passed = passed + 1;
else
    fprintf('  FAIL: radiusOfCurvature\n');
    failed = failed + 1;
end

% testGouyPhase
gouy = PhysicalConstants.gouyPhase(0.05, zr);
if (gouy > 0)
    fprintf('  PASS: gouyPhase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase\n');
    failed = failed + 1;
end

% test speed_of_light
c = PhysicalConstants.speed_of_light;
if (c > 1e8 && c < 3e8)
    fprintf('  PASS: speed_of_light\n');
    passed = passed + 1;
else
    fprintf('  FAIL: speed_of_light\n');
    failed = failed + 1;
end

% test planck
h = PhysicalConstants.planck;
if (h > 1e-34 && h < 1e-33)
    fprintf('  PASS: planck\n');
    passed = passed + 1;
else
    fprintf('  FAIL: planck\n');
    failed = failed + 1;
end

fprintf('\n=== PhysicalConstants: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
