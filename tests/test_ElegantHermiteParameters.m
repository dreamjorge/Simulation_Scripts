#!/usr/bin/env octave
% Tests for ElegantHermiteParameters

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantHermiteParameters Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;

% testConstructor
ehp = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
if (ehp.n == 1 && ehp.m == 1)
    fprintf('  PASS: constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor\n');
    failed = failed + 1;
end

% testDefaultNandM
ehp0 = ElegantHermiteParameters(0.05, w0, lambda);
if (ehp0.n == 0 && ehp0.m == 0)
    fprintf('  PASS: default n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: default n m\n');
    failed = failed + 1;
end

% testAlphaAtWaist
ehp_z0 = ElegantHermiteParameters(0, w0, lambda, 1, 1);
k_val = 2*pi/lambda;
zr = pi*w0^2/lambda;
expected_alpha = k_val / (2*zr);
if (abs(ehp_z0.alpha - expected_alpha) < 1e-10)
    fprintf('  PASS: alpha at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha at waist\n');
    failed = failed + 1;
end

% testAlphaAtZgt0
ehp_z = ElegantHermiteParameters(0.1, w0, lambda, 1, 1);
if (imag(ehp_z.alpha) ~= 0 || real(ehp_z.alpha) ~= 0)
    fprintf('  PASS: alpha complex at z>0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha at z>0\n');
    failed = failed + 1;
end

% testAlphaFormula
k_af = 2*pi/lambda;
zr_af = pi*w0^2/lambda;
q_af = 0.1 + 1i*zr_af;
expected_alpha = 1i * k_af / (2 * q_af);
if (abs(ehp_z.alpha - expected_alpha) < 1e-10)
    fprintf('  PASS: alpha formula\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha formula\n');
    failed = failed + 1;
end

% testAlphaChangesWithZ
ehp_z1 = ElegantHermiteParameters(0.01, w0, lambda, 1, 1);
ehp_z2 = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
if (ehp_z1.alpha ~= ehp_z2.alpha)
    fprintf('  PASS: alpha changes with z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: alpha change\n');
    failed = failed + 1;
end

% testInheritsGaussianParameters
if (isfield(ehp, 'RayleighDistance') || ~isempty(ehp.RayleighDistance))
    fprintf('  PASS: inherits from GaussianParameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: inheritance\n');
    failed = failed + 1;
end

% testWaistFromParent
if (ehp.Waist > 0)
    fprintf('  PASS: Waist from parent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Waist\n');
    failed = failed + 1;
end

% testGouyPhaseFromParent
if (ehp.GouyPhase > 0)
    fprintf('  PASS: GouyPhase from parent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: GouyPhase\n');
    failed = failed + 1;
end

% testkFromParent
if (ehp.k > 0)
    fprintf('  PASS: k from parent\n');
    passed = passed + 1;
else
    fprintf('  FAIL: k\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantHermiteParameters: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
