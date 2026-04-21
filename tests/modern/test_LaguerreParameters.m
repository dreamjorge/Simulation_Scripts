% Compatible with GNU Octave and MATLAB
% Tests for LaguerreParameters

repoRoot = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(fullfile(repoRoot, 'ParaxialBeams'));
addpath(fullfile(repoRoot, 'src', 'parameters'));
addpath(fullfile(repoRoot, 'src', 'computation'));

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
L_test = PolynomialUtils.associatedLaguerre(2, 1, x_test);
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

% testVectorZInput
z_vec = linspace(0, 0.1, 5);
lp_vec = LaguerreParameters(z_vec, w0, lambda, 1, 1);
if (numel(lp_vec.zCoordinate) == 5)
    fprintf('  PASS: vector z input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: vector z\n');
    failed = failed + 1;
end

% testStaticGetWaistAtWaist
wL_at_waist = LaguerreParameters.getWaist(0, w0, zr_s, 1, 1);
expected_lw = w0 * sqrt(2*1 + abs(1) + 1);
if (abs(wL_at_waist - expected_lw) / expected_lw < 1e-10)
    fprintf('  PASS: static getWaist at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getWaist at waist\n');
    failed = failed + 1;
end

% testInheritedProperties
lp_inherit = LaguerreParameters(z, w0, lambda, 1, 1);
if (lp_inherit.k > 0 && lp_inherit.RayleighDistance > 0 && lp_inherit.Waist > 0)
    fprintf('  PASS: inherited properties\n');
    passed = passed + 1;
else
    fprintf('  FAIL: inherited properties\n');
    failed = failed + 1;
end

% testHigherPLargerWaist
lp_p0 = LaguerreParameters(z, w0, lambda, 0, 0);
lp_p2 = LaguerreParameters(z, w0, lambda, 0, 2);
if (lp_p2.LaguerreWaist > lp_p0.LaguerreWaist)
    fprintf('  PASS: higher p larger waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher p\n');
    failed = failed + 1;
end

% testPhiPhaseWithNegativeL
lp_neg_phi = LaguerreParameters(z, w0, lambda, -2, 1);
expected_phi_neg = (abs(-2) + 2*1) * atan(z/zr);
if (abs(lp_neg_phi.PhiPhase - expected_phi_neg) < 1e-10)
    fprintf('  PASS: PhiPhase with negative l\n');
    passed = passed + 1;
else
    fprintf('  FAIL: PhiPhase negative l\n');
    failed = failed + 1;
end

% testLaguerreWaistFormula
lp_formula = LaguerreParameters(z, w0, lambda, 2, 3);
expected_formula = lp_formula.Waist * sqrt(2*3 + abs(2) + 1);
if (abs(lp_formula.LaguerreWaist - expected_formula) / expected_formula < 1e-10)
    fprintf('  PASS: LaguerreWaist formula\n');
    passed = passed + 1;
else
    fprintf('  FAIL: LaguerreWaist formula\n');
    failed = failed + 1;
end

% testAssociatedLaguerrePolynomialL0
L_l0 = PolynomialUtils.associatedLaguerre(1, 0, [0 1 2]);
if (all(isfinite(L_l0)))
    fprintf('  PASS: Laguerre polynomial L=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: L=0 polynomial\n');
    failed = failed + 1;
end

% testCombinedLandP
lp_comb = LaguerreParameters(z, w0, lambda, 3, 2);
expected_comb = lp_comb.Waist * sqrt(2*2 + abs(3) + 1);
if (abs(lp_comb.LaguerreWaist - expected_comb) / expected_comb < 1e-10)
    fprintf('  PASS: combined l p waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: combined l p\n');
    failed = failed + 1;
end

% --- Dynamic API tests (Phase 2) ---

% testDynamicPhiPhase
lp_dyn = LaguerreParameters(0, w0, lambda, 2, 1);
z2 = 0.1;
zr_dyn = lp_dyn.RayleighDistance;
phi_dyn = lp_dyn.phiPhase(z2);
phi_expected = (abs(2) + 2*1) * atan(z2 / zr_dyn);
if (abs(phi_dyn - phi_expected) < 1e-10)
    fprintf('  PASS: dynamic phiPhase(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dynamic phiPhase(z)\n');
    failed = failed + 1;
end

% testDynamicPhiPhaseAtZero
phi_zero = lp_dyn.phiPhase(0);
if (abs(phi_zero) < 1e-10)
    fprintf('  PASS: dynamic phiPhase(0) = 0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dynamic phiPhase(0)\n');
    failed = failed + 1;
end

% testDynamicLaguerreWaist
w_dyn2 = lp_dyn.laguerreWaist(z2);
w_g = lp_dyn.waist(z2);
w_expected2 = w_g * sqrt(2*1 + abs(2) + 1);
if (abs(w_dyn2 - w_expected2) / w_expected2 < 1e-10)
    fprintf('  PASS: dynamic laguerreWaist(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dynamic laguerreWaist(z)\n');
    failed = failed + 1;
end

% testDynamicPhiPhaseNegativeL
lp_neg = LaguerreParameters(0, w0, lambda, -2, 1);
phi_neg = lp_neg.phiPhase(z2);
phi_neg_exp = (abs(-2) + 2*1) * atan(z2 / lp_neg.RayleighDistance);
if (abs(phi_neg - phi_neg_exp) < 1e-10)
    fprintf('  PASS: dynamic phiPhase negative l\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dynamic phiPhase negative l\n');
    failed = failed + 1;
end

fprintf('\n=== LaguerreParameters: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
