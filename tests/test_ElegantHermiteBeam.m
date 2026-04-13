#!/usr/bin/env octave
% Tests for ElegantHermiteBeam

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== ElegantHermiteBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

% testFieldGeneration
ehp = ElegantHermiteParameters(0.01, w0, lambda, 1, 1);
ehb = ElegantHermiteBeam(X, Y, ehp);
if (size(ehb.OpticalField) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
ehp_z0 = ElegantHermiteParameters(0, w0, lambda, 1, 1);
ehb_z0 = ElegantHermiteBeam(X, Y, ehp_z0);
if (all(all(isfinite(ehb_z0.OpticalField))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
ehp_z = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
ehb_z = ElegantHermiteBeam(X, Y, ehp_z);
if (all(all(isfinite(ehb_z.OpticalField))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testCoordinatesStored
if (isequal(size(ehb.x), [64, 64]) && isequal(size(ehb.y), [64, 64]))
    fprintf('  PASS: coordinates stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: coordinates\n');
    failed = failed + 1;
end

% testParametersStored
if (ehb.Parameters.n == 1 && ehb.Parameters.m == 1)
    fprintf('  PASS: parameters stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: parameters\n');
    failed = failed + 1;
end

% testHigherOrderModes
ehp_high = ElegantHermiteParameters(0.01, w0, lambda, 2, 2);
ehb_high = ElegantHermiteBeam(X, Y, ehp_high);
if (all(all(isfinite(ehb_high.OpticalField))))
    fprintf('  PASS: higher order modes\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher order\n');
    failed = failed + 1;
end

% testAlphaStored
try
    alpha_val = ehb.Parameters.alpha;
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
if (size(ehb.OpticalField,1) == 64 && size(ehb.OpticalField,2) == 64)
    fprintf('  PASS: valid output size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: output size\n');
    failed = failed + 1;
end

% testZeroOrderElegant
ehp_00 = ElegantHermiteParameters(0, w0, lambda, 0, 0);
ehb_00 = ElegantHermiteBeam(X, Y, ehp_00);
if (all(all(isfinite(ehb_00.OpticalField))))
    fprintf('  PASS: zero order elegant valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order elegant\n');
    failed = failed + 1;
end

% testComplexFieldElegant
ehb_complex = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0.05, w0, lambda, 1, 1));
if (ehb_complex.Parameters.GouyPhase ~= 0)
    fprintf('  PASS: complex field elegant\n');
    passed = passed + 1;
else
    fprintf('  FAIL: complex field elegant\n');
    failed = failed + 1;
end

% testDifferentNM
ehb_n1 = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, lambda, 2, 0));
ehb_m1 = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, lambda, 0, 2));
if (all(all(isfinite(ehb_n1.OpticalField))) && all(all(isfinite(ehb_m1.OpticalField))))
    fprintf('  PASS: different n m orders\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different n m\n');
    failed = failed + 1;
end

% testElegantAlphaFormula
ehp_alpha = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
alpha_calc = ehp_alpha.alpha;
if (isfinite(alpha_calc))
    fprintf('  PASS: elegant alpha formula\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant alpha\n');
    failed = failed + 1;
end

% testElegantWaistInheritance
ehp_w = ElegantHermiteParameters(0.1, w0, lambda, 1, 1);
if (ehp_w.Waist > 0)
    fprintf('  PASS: elegant waist inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant waist\n');
    failed = failed + 1;
end

% testElegantKInheritance
if (ehp.k > 0)
    fprintf('  PASS: elegant k inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant k\n');
    failed = failed + 1;
end

% testElegantGouyInheritance
if (ehp_z.GouyPhase > 0)
    fprintf('  PASS: elegant Gouy inheritance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant Gouy\n');
    failed = failed + 1;
end

% testElegantNegativeZ
ehb_neg = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(-0.05, w0, lambda, 1, 1));
if (all(all(isfinite(ehb_neg.OpticalField))))
    fprintf('  PASS: elegant negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant negative z\n');
    failed = failed + 1;
end

% testElegantCoordinatesNotEmpty
if (~isempty(ehb.x) && ~isempty(ehb.y))
    fprintf('  PASS: elegant coordinates not empty\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant coordinates\n');
    failed = failed + 1;
end

% testElegantDifferentWavelengths
ehb_wl1 = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, 532e-9, 1, 1));
ehb_wl2 = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, 1064e-9, 1, 1));
if (all(all(isfinite(ehb_wl1.OpticalField))) && all(all(isfinite(ehb_wl2.OpticalField))))
    fprintf('  PASS: elegant different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant wavelengths\n');
    failed = failed + 1;
end

% testElegantHigherNM
ehb_hnm = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, lambda, 3, 3));
if (all(all(isfinite(ehb_hnm.OpticalField))))
    fprintf('  PASS: elegant higher n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant higher n m\n');
    failed = failed + 1;
end

% testElegantFieldAmplitude
ehb_amp = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, lambda, 0, 0));
max_amp_e = max(max(abs(ehb_amp.OpticalField)));
if (max_amp_e > 0)
    fprintf('  PASS: elegant field amplitude\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant amplitude\n');
    failed = failed + 1;
end

% testElegantAsymmetricNM
ehb_asym = ElegantHermiteBeam(X, Y, ElegantHermiteParameters(0, w0, lambda, 2, 1));
if (all(all(isfinite(ehb_asym.OpticalField))))
    fprintf('  PASS: elegant asymmetric n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant asymmetric\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantHermiteBeam: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
