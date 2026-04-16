% Compatible with GNU Octave and MATLAB
% Edge case tests for HankelLaguerre
%
% Edge cases covered:
%   - r = 0 (origin in polar coordinates)
%   - theta = 0 (azimuthal angle)
%   - z = 0 (beam waist plane)
%   - l = 0 (azimuthal order zero)
%   - p = 0 (radial order zero)
%   - l and p both zero (fundamental mode)
%
% RISKS IDENTIFIED:
%   1. r = 0: amp = (sqrt(2)*r/w).^|l| = 0 for l>0, but multiplied by Lpl which is 1 at origin
%   2. r = 0: theta undefined but never used when r=0 because amp=0
%   3. z = 0: Rc=Inf handled by: Rc(z==0)=Inf; phase_curv(isinf(Rc))=0;
%   4. l = 0: exp(i*l*theta) = 1, no angular dependence
%   5. H^(1) and H^(2) Hankel types: superposition of LB + i*XLG or LB - i*XLG

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== HankelLaguerre Edge Case Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 1e-3;
lambda = 632e-9;
grid = GridUtils(64, 64, 4e-3, 4e-3);
[X, Y] = grid.create2DGrid();
[R, TH] = cart2pol(X, Y);

% test_z0_Type1
% RISK: z=0 with HankelType=1 (H^1)
HL1 = HankelLaguerre(w0, lambda, 1, 0, 1);
field1 = HL1.opticalField(R, TH, 0);
if (all(all(isfinite(field1))))
    fprintf('  PASS: HankelType=1 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=1 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_Type2
% RISK: z=0 with HankelType=2 (H^2)
HL2 = HankelLaguerre(w0, lambda, 1, 0, 2);
field2 = HL2.opticalField(R, TH, 0);
if (all(all(isfinite(field2))))
    fprintf('  PASS: HankelType=2 at z=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType=2 at z=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_r0_l1
% RISK: r=0 with l=1, amp=(sqrt(2)*r/w)^1=0 so field should be 0 at origin
% But the amplitude is 0 at origin for l>0 modes (on-axis intensity null)
HL_l1 = HankelLaguerre(w0, lambda, 1, 0, 1);
field_l1 = HL_l1.opticalField(R, TH, 0.01);
center_val = field_l1(33,33);  % center of 64x64 grid
if (isfinite(center_val) && abs(center_val) < 0.1)
    fprintf('  PASS: l=1 mode has null at origin (expected ~0)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: l=1 mode should have null at origin, got %g\n', abs(center_val));
    failed = failed + 1;
end

% test_r0_l0
% RISK: r=0 with l=0, amp=(sqrt(2)*r/w)^0=1, field finite at origin
HL_l0 = HankelLaguerre(w0, lambda, 0, 0, 1);
field_l0 = HL_l0.opticalField(R, TH, 0.01);
center_val_l0 = field_l0(33,33);
if (isfinite(center_val_l0) && abs(center_val_l0) > 0.5)
    fprintf('  PASS: l=0 mode has maximum at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: l=0 mode should have maximum at origin, got %g\n', abs(center_val_l0));
    failed = failed + 1;
end

% test_theta0
% RISK: theta=0 should not cause issues (exp(i*l*0)=1)
HL_theta0 = HankelLaguerre(w0, lambda, 1, 0, 1);
field_theta0 = HL_theta0.opticalField(R, TH, 0);
if (all(all(isfinite(field_theta0))))
    fprintf('  PASS: theta=0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: theta=0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_z0_l0p0
% RISK: z=0 with l=0, p=0 (fundamental mode)
HL_fund = HankelLaguerre(w0, lambda, 0, 0, 1);
field_fund = HL_fund.opticalField(R, TH, 0);
if (all(all(isfinite(field_fund))) && abs(field_fund(33,33)) > 0.5)
    fprintf('  PASS: fundamental mode (l=0,p=0) finite at origin z=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fundamental mode issue at origin z=0\n');
    failed = failed + 1;
end

% test_z0_SameAs_eps
% Verify z=0 and z=eps produce approximately same result
field_fund_eps = HL_fund.opticalField(R, TH, eps);
if (max(max(abs(field_fund - field_fund_eps))) < 1e-6)
    fprintf('  PASS: z=0 and z=eps produce same field for fundamental\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=0 and z=eps differ for fundamental mode\n');
    failed = failed + 1;
end

% test_pNonZero
% RISK: p>0 adds radial nodes, should still be finite
HL_p1 = HankelLaguerre(w0, lambda, 0, 1, 1);
field_p1 = HL_p1.opticalField(R, TH, 0);
if (all(all(isfinite(field_p1))))
    fprintf('  PASS: p=1 (one radial node) produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: p=1 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_zVeryLarge
% RISK: Very large z may cause overflow in z/zr ratio
field_zlarge = HL1.opticalField(R, TH, 1e10*w0);
if (all(all(isfinite(field_zlarge))))
    fprintf('  PASS: z=1e10*w0 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=1e10*w0 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_zNegative
% z negative (before waist) should be valid
field_zneg = HL1.opticalField(R, TH, -0.01);
if (all(all(isfinite(field_zneg))))
    fprintf('  PASS: z=-0.01 produces finite field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z=-0.01 produces NaN/Inf\n');
    failed = failed + 1;
end

% test_H1_vs_H2_Difference
% H^1 and H^2 should be different (complex conjugate relationship)
diff_H1_H2 = max(max(abs(field1 - field2)));
if (diff_H1_H2 > 1e-10)
    fprintf('  PASS: H^1 and H^2 produce different fields (as expected)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: H^1 and H^2 produce same field (unexpected)\n');
    failed = failed + 1;
end

fprintf('\n=== HankelLaguerre Edge Case Summary ===\n');
fprintf('Passed: %d\n', passed);
fprintf('Failed: %d\n', failed);

if failed > 0
    fprintf('ESTADO: FALLO\n');
else
    fprintf('ESTADO: ÉXITO\n');
end