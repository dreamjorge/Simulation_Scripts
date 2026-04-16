% Compatible with GNU Octave and MATLAB
% Edge case tests for HankelHermite
%
% Edge cases covered:
%   - z = 0 (beam waist plane)
%   - Origin (0,0) in Cartesian grid (r=0)
%   - r=0 in Hankel functions (Bessel singular behavior)
%   - theta=0 (azimuthal angle)
%   - Extreme orders (n=0, m=0, n>0, m>0)
%
% RISKS IDENTIFIED:
%   1. r=0 at origin: Hankel functions are well-behaved but Bessel J_n(0)=0 for n>0
%   2. z=0: Rc=Inf handled by: Rc(z==0)=Inf; phase_curv(isinf(Rc))=0;
%   3. theta undefined at r=0 but multiplied by r^|l| so mathematically OK
%   4. NHG (NHx, NHy) solutions must not diverge at r=0

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== HankelHermite Edge Case Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 1e-3;
lambda = 632e-9;
grid = GridUtils(64, 64, 4e-3, 4e-3);
[X, Y] = grid.create2DGrid();

% test_z0_Type11
% RISK: z=0 with HankelType=11
HH11 = HankelHermite(w0, lambda, 1, 1, 11);
field11 = HH11.opticalField(X, Y, 0);
if (all(all(isfinite(field11))))
    fprintf('  PASS: HankelType=11 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=11 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_Type12
HH12 = HankelHermite(w0, lambda, 1, 1, 12);
field12 = HH12.opticalField(X, Y, 0);
if (all(all(isfinite(field12))))
    fprintf('  PASS: HankelType=12 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=12 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_Type21
HH21 = HankelHermite(w0, lambda, 1, 1, 21);
field21 = HH21.opticalField(X, Y, 0);
if (all(all(isfinite(field21))))
    fprintf('  PASS: HankelType=21 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=21 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_Type22
HH22 = HankelHermite(w0, lambda, 1, 1, 22);
field22 = HH22.opticalField(X, Y, 0);
if (all(all(isfinite(field22))))
    fprintf('  PASS: HankelType=22 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=22 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_Origin_r0
% RISK: r=0 at origin (center of grid) - Hankel functions well-behaved
r = sqrt(X.^2 + Y.^2);
origin_idx = 33;  % center of 64x64 grid
field_origin = HH11.opticalField(X, Y, 0.01);
origin_val = field_origin(origin_idx, origin_idx);
if (isfinite(origin_val))
    fprintf('  PASS: field finite at origin r=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at origin r=0 is NaN/Inf\n');
    failed = failed + 1;
end

% test_ZeroOrder_n0m0
% RISK: Zero-order (n=0, m=0) should reduce to Gaussian-like behavior
HH00 = HankelHermite(w0, lambda, 0, 0, 11);
field00 = HH00.opticalField(X, Y, 0);
if (all(all(isfinite(field00))) && abs(field00(origin_idx,origin_idx)) > 0.5)
    fprintf('  PASS: n=0,m=0 produces finite significant field at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: n=0,m=0 field issue\n');
    failed = failed + 1;
end

% test_HighOrder_n5m5
% RISK: High orders may cause numerical issues in Hermite polynomials
HH55 = HankelHermite(w0, lambda, 5, 5, 11);
field55 = HH55.opticalField(X, Y, 0);
if (all(all(isfinite(field55))))
    fprintf('  PASS: n=5,m=5 high order produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: n=5,m=5 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_zVeryLarge
% RISK: Very large z may cause overflow in z/zr ratio
field_zlarge = HH11.opticalField(X, Y, 1e10*w0);
if (all(all(isfinite(field_zlarge))))
    fprintf('  PASS: z=1e10*w0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=1e10*w0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_zNegative
% z negative (before waist) should be valid
field_zneg = HH11.opticalField(X, Y, -0.01);
if (all(all(isfinite(field_zneg))))
    fprintf('  PASS: z=-0.01 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=-0.01 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_NHG_atOrigin
% RISK: NHG (normalized Hermite) should be finite at origin
% The Hankel combination uses (Hx +/- i*NHx)*(Hy +/- i*NHy)
% If NHG diverges at origin, Hankel fields would be infinite
[Hx0, NHx0] = HermiteParameters.getHermiteSolutions(0, 0);
[Hx1, NHx1] = HermiteParameters.getHermiteSolutions(1, 0);
if (isfinite(NHx0) && isfinite(NHx1))
    fprintf('  PASS: NHG solutions finite at origin for n=0,1\n');
    passed = passed + 1;
else
    fprintf('  FAIL: NHG diverges at origin\n');
    failed = failed + 1;
end

fprintf('\n=== HankelHermite Edge Case Summary ===\n');
fprintf('Passed: %d\n', passed);
fprintf('Failed: %d\n', failed);

if failed > 0
    fprintf('ESTADO: FALLO\n');
else
    fprintf('ESTADO: ÉXITO\n');
end