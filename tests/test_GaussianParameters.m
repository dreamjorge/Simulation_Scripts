#!/usr/bin/env octave
% Tests for GaussianParameters

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== GaussianParameters Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
z = 0.05;

% testRayleighDistance
params = GaussianParameters(z, w0, lambda);
zr_expected = pi * w0^2 / lambda;
if (abs(params.RayleighDistance - zr_expected) / zr_expected < 1e-5)
    fprintf('  PASS: RayleighDistance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: RayleighDistance\n');
    failed = failed + 1;
end

% testWaistAtOrigin
if (params.Waist >= w0)
    fprintf('  PASS: waist at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist at origin\n');
    failed = failed + 1;
end

% testWaveNumber
expected_k = 2*pi / lambda;
if (abs(params.k - expected_k) / expected_k < 1e-5)
    fprintf('  PASS: waveNumber k\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber k\n');
    failed = failed + 1;
end

% testVectorZ
z_vec = linspace(0, 0.1, 10);
params_vec = GaussianParameters(z_vec, w0, lambda);
if (numel(params_vec.zCoordinate) == 10)
    fprintf('  PASS: vector z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: vector z\n');
    failed = failed + 1;
end

% testWaistIncreasesWithZ
params1 = GaussianParameters(0, w0, lambda);
params2 = GaussianParameters(0.1, w0, lambda);
if (params2.Waist > params1.Waist)
    fprintf('  PASS: waist increases with z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist increases with z\n');
    failed = failed + 1;
end

% testGouyPhase
if (params.GouyPhase > 0)
    fprintf('  PASS: GouyPhase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: GouyPhase\n');
    failed = failed + 1;
end

% testRadiusAtWaist
params0 = GaussianParameters(0, w0, lambda);
if (isinf(params0.Radius))
    fprintf('  PASS: Radius at waist is Inf\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Radius at waist\n');
    failed = failed + 1;
end

% testDivergenceAngle
if (params.DivergenceAngle > 0)
    fprintf('  PASS: DivergenceAngle\n');
    passed = passed + 1;
else
    fprintf('  FAIL: DivergenceAngle\n');
    failed = failed + 1;
end

% testToString
str = params.toString();
if (~isempty(str))
    fprintf('  PASS: toString\n');
    passed = passed + 1;
else
    fprintf('  FAIL: toString\n');
    failed = failed + 1;
end

% testIsEqual
p1 = GaussianParameters(0, w0, lambda);
p2 = GaussianParameters(0, w0, lambda);
if (p1.isEqual(p2))
    fprintf('  PASS: isEqual\n');
    passed = passed + 1;
else
    fprintf('  FAIL: isEqual\n');
    failed = failed + 1;
end

% testStaticRayleighDistance
zr_s = GaussianParameters.rayleighDistance(w0, lambda);
if (abs(zr_s - zr_expected) / zr_expected < 1e-5)
    fprintf('  PASS: static rayleighDistance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: static rayleighDistance\n');
    failed = failed + 1;
end

% testStaticGetWaist
w_s = GaussianParameters.getWaist(0.05, w0, zr_expected);
if (w_s > 0)
    fprintf('  PASS: static getWaist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: static getWaist\n');
    failed = failed + 1;
end

fprintf('\n=== GaussianParameters: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
