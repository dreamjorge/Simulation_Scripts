#!/usr/bin/env octave
% Tests for ElegantLaguerreBeam (Phase 3 API: ElegantLaguerreBeam(w0, lambda, l, p))

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantLaguerreBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

elb = ElegantLaguerreBeam(w0, lambda, 1, 0);

% testFieldGeneration
field = elb.opticalField(X, Y, 0.01);
if (size(field) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
field_z0 = elb.opticalField(X, Y, 0);
if (all(all(isfinite(field_z0))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
field_z = elb.opticalField(X, Y, 0.05);
if (all(all(isfinite(field_z))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testModeIndicesStored
if (elb.l == 1 && elb.p == 0)
    fprintf('  PASS: mode indices stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: mode indices\n');
    failed = failed + 1;
end

% testInitialWaistStored
if (elb.InitialWaist == w0)
    fprintf('  PASS: InitialWaist stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: InitialWaist\n');
    failed = failed + 1;
end

% testGetParameters
params = elb.getParameters(0.05);
if (abs(params.zCoordinate - 0.05) < 1e-15 && params.InitialWaist == w0)
    fprintf('  PASS: getParameters(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getParameters(z)\n');
    failed = failed + 1;
end

% testHigherOrderModes
elb_high = ElegantLaguerreBeam(w0, lambda, 2, 1);
field_high = elb_high.opticalField(X, Y, 0.01);
if (all(all(isfinite(field_high))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order\n');
    failed = failed + 1;
end

% testValidOutputSize
if (size(field_z0, 1) == 64 && size(field_z0, 2) == 64)
    fprintf('  PASS: valid output size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: output size\n');
    failed = failed + 1;
end

% testDifferentLandP
elb_10 = ElegantLaguerreBeam(w0, lambda, 1, 0);
elb_01 = ElegantLaguerreBeam(w0, lambda, 0, 1);
f_10 = elb_10.opticalField(X, Y, 0.01);
f_01 = elb_01.opticalField(X, Y, 0.01);
if (size(f_10) == size(f_01))
    fprintf('  PASS: different l p valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different l p\n');
    failed = failed + 1;
end

% testZeroOrderElegantLaguerre
elb_00 = ElegantLaguerreBeam(w0, lambda, 0, 0);
field_00 = elb_00.opticalField(X, Y, 0);
if (all(all(isfinite(field_00))))
    fprintf('  PASS: zero order elegant LG valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order elegant LG\n');
    failed = failed + 1;
end

% testGouyPhaseIncluded
params_z = elb.getParameters(0.05);
if (params_z.GouyPhase > 0)
    fprintf('  PASS: Gouy phase included\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Gouy phase\n');
    failed = failed + 1;
end

% testHigherPOrder
elb_p2 = ElegantLaguerreBeam(w0, lambda, 0, 2);
field_p2 = elb_p2.opticalField(X, Y, 0);
if (all(all(isfinite(field_p2))))
    fprintf('  PASS: higher p order elegant\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher p order elegant\n');
    failed = failed + 1;
end

% testNegativeLElegant
elb_neg = ElegantLaguerreBeam(w0, lambda, -1, 0);
field_neg = elb_neg.opticalField(X, Y, 0);
if (all(all(isfinite(field_neg))))
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

% testWaistFromParameters
params_w = elb.getParameters(0.1);
if (params_w.Waist > 0)
    fprintf('  PASS: waist from parameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist from parameters\n');
    failed = failed + 1;
end

% testElegantLGNegativeZ
field_neg_z = elb.opticalField(X, Y, -0.05);
if (all(all(isfinite(field_neg_z))))
    fprintf('  PASS: elegant LG negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG negative z\n');
    failed = failed + 1;
end

% testElegantLGDifferentWavelengths
elb_wl1 = ElegantLaguerreBeam(w0, 532e-9, 1, 0);
elb_wl2 = ElegantLaguerreBeam(w0, 1064e-9, 1, 0);
f_wl1 = elb_wl1.opticalField(X, Y, 0);
f_wl2 = elb_wl2.opticalField(X, Y, 0);
if (all(all(isfinite(f_wl1))) && all(all(isfinite(f_wl2))))
    fprintf('  PASS: elegant LG different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG wavelengths\n');
    failed = failed + 1;
end

% testElegantLGCombinedLP
elb_comb = ElegantLaguerreBeam(w0, lambda, 2, 3);
f_comb = elb_comb.opticalField(X, Y, 0);
if (all(all(isfinite(f_comb))))
    fprintf('  PASS: elegant LG combined l p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG combined\n');
    failed = failed + 1;
end

% testElegantLGFieldAmplitude
field_amp = elb_00.opticalField(X, Y, 0);
max_amp_el = max(max(abs(field_amp)));
if (max_amp_el > 0)
    fprintf('  PASS: elegant LG field amplitude\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG amplitude\n');
    failed = failed + 1;
end

% testElegantLGNegativeLWithP
elb_negp = ElegantLaguerreBeam(w0, lambda, -2, 1);
f_negp = elb_negp.opticalField(X, Y, 0);
if (all(all(isfinite(f_negp))))
    fprintf('  PASS: elegant LG negative l with p\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant LG negative l with p\n');
    failed = failed + 1;
end

% testFieldAtDifferentZ
f_z1 = elb.opticalField(X, Y, 0.01);
f_z2 = elb.opticalField(X, Y, 0.1);
if (all(all(isfinite(f_z1))) && all(all(isfinite(f_z2))))
    fprintf('  PASS: field at different z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at different z\n');
    failed = failed + 1;
end

% testBeamName
if (strcmp(elb.beamName(), 'elegant_laguerre_1_0'))
    fprintf('  PASS: beamName\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beamName\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantLaguerreBeam: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
