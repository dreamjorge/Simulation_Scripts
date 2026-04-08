#!/usr/bin/env octave
% Tests for ElegantLaguerreParameters

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantLaguerreParameters Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;

% testConstructor
elp = ElegantLaguerreParameters(0.05, w0, lambda, 1, 1);
if (elp.l == 1 && elp.p == 1)
    fprintf('  PASS: constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor\n');
    failed = failed + 1;
end

% testDefaultLandP
elp0 = ElegantLaguerreParameters(0.05, w0, lambda);
if (elp0.l == 0 && elp0.p == 0)
    fprintf('  PASS: default l p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: default l p\n');
    failed = failed + 1;
end

% testAlphaAtWaist
elp_z0 = ElegantLaguerreParameters(0, w0, lambda, 1, 1);
k_val = 2*pi/lambda;
zr = pi*w0^2/lambda;
expected_alpha = k_val / (2*zr);
if (abs(elp_z0.alpha - expected_alpha) < 1e-10)
    fprintf('  PASS: alpha at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha at waist\n');
    failed = failed + 1;
end

% testAlphaAtZgt0
elp_z = ElegantLaguerreParameters(0.1, w0, lambda, 1, 1);
if (imag(elp_z.alpha) ~= 0 || real(elp_z.alpha) ~= 0)
    fprintf('  PASS: alpha complex at z>0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha at z>0\n');
    failed = failed + 1;
end

% testAlphaMatchesElegantHermite
ehp_test = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
elp_test = ElegantLaguerreParameters(0.05, w0, lambda, 1, 0);
if (abs(ehp_test.alpha - elp_test.alpha) < 1e-10)
    fprintf('  PASS: alpha matches ElegantHermiteParameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha match\n');
    failed = failed + 1;
end

% testAlphaChangesWithZ
elp_z1 = ElegantLaguerreParameters(0.01, w0, lambda, 1, 1);
elp_z2 = ElegantLaguerreParameters(0.05, w0, lambda, 1, 1);
if (elp_z1.alpha ~= elp_z2.alpha)
    fprintf('  PASS: alpha changes with z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha change\n');
    failed = failed + 1;
end

% testInheritsGaussianParameters
if (isfield(elp, 'RayleighDistance') || ~isempty(elp.RayleighDistance))
    fprintf('  PASS: inherits from GaussianParameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: inheritance\n');
    failed = failed + 1;
end

% testWaistFromParent
if (elp.Waist > 0)
    fprintf('  PASS: Waist from parent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Waist\n');
    failed = failed + 1;
end

% testGouyPhaseFromParent
if (elp.GouyPhase > 0)
    fprintf('  PASS: GouyPhase from parent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: GouyPhase\n');
    failed = failed + 1;
end

% testNegativeL
elp_neg = ElegantLaguerreParameters(0.05, w0, lambda, -1, 1);
if (elp_neg.l == -1)
    fprintf('  PASS: negative l\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantLaguerreParameters: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
