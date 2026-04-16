% Compatible with GNU Octave and MATLAB
% Edge case tests for GaussianBeam
%
% Edge cases covered:
%   - z = 0 (beam waist plane)
%   - Origin (0,0) in Cartesian coordinates
%   - Extreme parameters (w0 near eps, w0 very large, z very large)
%
% RISKS IDENTIFIED:
%   1. z = 0 causes Rc = Inf in curvature term, handled by: Rc(z==0) = Inf; phase_curv(isinf(Rc)) = 0;
%   2. Origin r=0 causes no division by zero in Hankel functions (not used by Gaussian)
%   3. w0 near eps may cause numerical underflow
%   4. z very large may cause overflow in z/zr ratio

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== GaussianBeam Edge Case Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

beam = GaussianBeam(w0, lambda);

% test_z0_NoNaN
% RISK: z=0 could cause division by zero in Rc calculation
field_z0 = beam.opticalField(X, Y, 0);
if (all(all(isfinite(field_z0))))
    fprintf('  PASS: z=0 produces finite field (no NaN/Inf)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_AmplitudeAtCenter
% At waist (z=0), amplitude should be maximum at origin
field_center = beam.opticalField(zeros(64,64), zeros(64,64), 0);
if (abs(abs(field_center(33,33)) - 1) < 1e-10)
    fprintf('  PASS: amplitude at origin z=0 is normalized to 1\n');
    passed = passed + 1;
else
    fprintf('  FAIL: amplitude at origin z=0 not normalized\n');
    failed = failed + 1;
end

% test_z0_PhaseAtOrigin
% At waist, Gouy phase should be zero
phase_center = angle(field_center(33,33));
if (abs(phase_center) < 1e-10)
    fprintf('  PASS: phase at origin z=0 is zero (Gouy=0 at waist)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase at origin z=0 not zero: %g\n', phase_center);
    failed = failed + 1;
end

% test_z0_SameAs_eps
% Verify z=0 and z=eps produce approximately same result
field_eps = beam.opticalField(X, Y, eps);
if (max(max(abs(field_z0 - field_eps))) < 1e-8)
    fprintf('  PASS: z=0 and z=eps produce same field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=0 and z=eps differ\n');
    failed = failed + 1;
end

% test_Origin_ExactZero
% RISK: Origin (0,0) exact zeros may cause issues in r^2/w^2
field_origin = beam.opticalField(X, Y, 0.01);
center_val = field_origin(33,33);
if (isfinite(center_val) && abs(center_val) > 0.5)
    fprintf('  PASS: field at origin is finite and significant\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at origin issue: finite=%d, value=%g\n', isfinite(center_val), abs(center_val));
    failed = failed + 1;
end

% test_w0_VerySmall
% RISK: w0 near eps may cause numerical underflow in amplitude
beam_small = GaussianBeam(eps, lambda);
field_small = beam_small.opticalField(X, Y, 0);
if (all(all(isfinite(field_small))))
    fprintf('  PASS: w0=eps produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: w0=eps produces NaN/Inf\n');
    failed = failed + 1;
end

% test_w0_VeryLarge
% RISK: w0 very large may cause numerical overflow
beam_large = GaussianBeam(1e10*w0, lambda);
field_large = beam_large.opticalField(X, Y, 0);
if (all(all(isfinite(field_large))))
    fprintf('  PASS: w0=1e10*w0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: w0=1e10*w0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z_VeryLarge
% RISK: z very large may cause overflow in z/zr ratio
beam_zlarge = GaussianBeam(w0, lambda);
field_zlarge = beam_zlarge.opticalField(X, Y, 1e10*w0);
if (all(all(isfinite(field_zlarge))))
    fprintf('  PASS: z=1e10*w0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=1e10*w0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z_Negative
% z negative should be valid (before waist)
field_zneg = beam.opticalField(X, Y, -0.01);
if (all(all(isfinite(field_zneg))))
    fprintf('  PASS: z=-0.01 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=-0.01 produces NaN/Inf\n');
    failed = failed + 1;
end

fprintf('\n=== GaussianBeam Edge Case Summary ===\n');
fprintf('Passed: %d\n', passed);
fprintf('Failed: %d\n', failed);

if failed > 0
    fprintf('ESTADO: FALLO\n');
else
    fprintf('ESTADO: ÉXITO\n');
end