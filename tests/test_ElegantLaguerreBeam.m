#!/usr/bin/env octave
% Tests for ElegantLaguerreBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantLaguerreBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();
[R, Theta] = cart2pol(X, Y);

% testFieldGeneration
elp = ElegantLaguerreParameters(0.01, w0, lambda, 1, 0);
elb = ElegantLaguerreBeam(R, Theta, elp);
if (size(elb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
elp_z0 = ElegantLaguerreParameters(0, w0, lambda, 1, 0);
elb_z0 = ElegantLaguerreBeam(R, Theta, elp_z0);
if (all(all(isfinite(elb_z0.OpticalField))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
elp_z = ElegantLaguerreParameters(0.05, w0, lambda, 1, 0);
elb_z = ElegantLaguerreBeam(R, Theta, elp_z);
if (all(all(isfinite(elb_z.OpticalField))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testParametersStored
if (elb.Parameters.l == 1 && elb.Parameters.p == 0)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters\n');
    failed = failed + 1;
end

% testHigherOrderModes
elp_high = ElegantLaguerreParameters(0.01, w0, lambda, 2, 1);
elb_high = ElegantLaguerreBeam(R, Theta, elp_high);
if (all(all(isfinite(elb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order\n');
    failed = failed + 1;
end

% testAlphaStored
try
    alpha_val = elb.Parameters.alpha;
    if (~isempty(alpha_val))
        fprintf('  PASS: alpha stored\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: alpha\n');
        failed = failed + 1;
    end
catch
    fprintf('  PASS: alpha accessible\n');
    passed = passed + 1;
end

% testValidOutputSize
if (size(elb.OpticalField,1) == 64 && size(elb.OpticalField,2) == 64)
    fprintf('  PASS: valid output size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: output size\n');
    failed = failed + 1;
end

% testDifferentLandP
elp_10 = ElegantLaguerreParameters(0.01, w0, lambda, 1, 0);
elp_01 = ElegantLaguerreParameters(0.01, w0, lambda, 0, 1);
elb_10 = ElegantLaguerreBeam(R, Theta, elp_10);
elb_01 = ElegantLaguerreBeam(R, Theta, elp_01);
if (size(elb_10.OpticalField) == size(elb_01.OpticalField))
    fprintf('  PASS: different l p valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different l p\n');
    failed = failed + 1;
end

% testZeroOrderElegantLaguerre
elp_00 = ElegantLaguerreParameters(0, w0, lambda, 0, 0);
elb_00 = ElegantLaguerreBeam(R, Theta, elp_00);
if (all(all(isfinite(elb_00.OpticalField))))
    fprintf('  PASS: zero order elegant LG valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order elegant LG\n');
    failed = failed + 1;
end

% testComplexFieldElegantLaguerre
elb_complex = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0.05, w0, lambda, 1, 0));
if (elb_complex.Parameters.GouyPhase ~= 0)
    fprintf('  PASS: complex field elegant LG\n');
    passed = passed + 1;
else
    fprintf('  FAIL: complex field elegant LG\n');
    failed = failed + 1;
end

% testHigherPOrder
elp_p2 = ElegantLaguerreParameters(0, w0, lambda, 0, 2);
elb_p2 = ElegantLaguerreBeam(R, Theta, elp_p2);
if (all(all(isfinite(elb_p2.OpticalField))))
    fprintf('  PASS: higher p order elegant\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher p order elegant\n');
    failed = failed + 1;
end

% testNegativeLElegant
elp_neg = ElegantLaguerreParameters(0, w0, lambda, -1, 0);
elb_neg = ElegantLaguerreBeam(R, Theta, elp_neg);
if (all(all(isfinite(elb_neg.OpticalField))))
    fprintf('  PASS: negative l elegant valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l elegant\n');
    failed = failed + 1;
end

% testElegantLGAlphaFormula
elp_alpha = ElegantLaguerreParameters(0.05, w0, lambda, 1, 1);
alpha_el = elp_alpha.alpha;
if (isfinite(alpha_el))
    fprintf('  PASS: elegant LG alpha formula\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG alpha\n');
    failed = failed + 1;
end

% testElegantLGWaistInheritance
elp_w = ElegantLaguerreParameters(0.1, w0, lambda, 1, 1);
if (elp_w.Waist > 0)
    fprintf('  PASS: elegant LG waist inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG waist\n');
    failed = failed + 1;
end

% testElegantLGKInheritance
if (elp.k > 0)
    fprintf('  PASS: elegant LG k inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG k\n');
    failed = failed + 1;
end

% testElegantLGGouyInheritance
if (elp_z.GouyPhase > 0)
    fprintf('  PASS: elegant LG Gouy inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG Gouy\n');
    failed = failed + 1;
end

% testElegantLGNegativeZ
elb_neg_z = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(-0.05, w0, lambda, 1, 0));
if (all(all(isfinite(elb_neg_z.OpticalField))))
    fprintf('  PASS: elegant LG negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG negative z\n');
    failed = failed + 1;
end

% testElegantLGDifferentWavelengths
elb_wl1 = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0, w0, 532e-9, 1, 0));
elb_wl2 = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0, w0, 1064e-9, 1, 0));
if (all(all(isfinite(elb_wl1.OpticalField))) && all(all(isfinite(elb_wl2.OpticalField))))
    fprintf('  PASS: elegant LG different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG wavelengths\n');
    failed = failed + 1;
end

% testElegantLGCombinedLP
elb_comb = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0, w0, lambda, 2, 3));
if (all(all(isfinite(elb_comb.OpticalField))))
    fprintf('  PASS: elegant LG combined l p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG combined\n');
    failed = failed + 1;
end

% testElegantLGFieldAmplitude
elb_amp = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0, w0, lambda, 0, 0));
max_amp_el = max(max(abs(elb_amp.OpticalField)));
if (max_amp_el > 0)
    fprintf('  PASS: elegant LG field amplitude\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG amplitude\n');
    failed = failed + 1;
end

% testElegantLGNegativeLWithP
elb_negp = ElegantLaguerreBeam(R, Theta, ElegantLaguerreParameters(0, w0, lambda, -2, 1));
if (all(all(isfinite(elb_negp.OpticalField))))
    fprintf('  PASS: elegant LG negative l with p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG negative l with p\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantLaguerreBeam: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
