% Compatible with GNU Octave and MATLAB
% Tests for ElegantHermiteBeam (Phase 3 API: ElegantHermiteBeam(w0, lambda, n, m))

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

fprintf('=== ElegantHermiteBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

ehb = ElegantHermiteBeam(w0, lambda, 1, 1);

% testFieldGeneration
field = ehb.opticalField(X, Y, 0.01);
if (isequal(size(field), [64, 64]))
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAtWaist
field_z0 = ehb.opticalField(X, Y, 0);
if (all(all(isfinite(field_z0))))
    fprintf('  PASS: at waist finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at waist\n');
    failed = failed + 1;
end

% testAtZgt0
field_z = ehb.opticalField(X, Y, 0.05);
if (all(all(isfinite(field_z))))
    fprintf('  PASS: at z>0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z>0\n');
    failed = failed + 1;
end

% testModeOrdersStored
if (ehb.n == 1 && ehb.m == 1)
    fprintf('  PASS: mode orders stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: mode orders\n');
    failed = failed + 1;
end

% testInitialWaistStored
if (ehb.InitialWaist == w0)
    fprintf('  PASS: InitialWaist stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: InitialWaist\n');
    failed = failed + 1;
end

% testGetParameters
params = ehb.getParameters(0.05);
if (abs(params.zCoordinate - 0.05) < 1e-15 && params.InitialWaist == w0)
    fprintf('  PASS: getParameters(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getParameters(z)\n');
    failed = failed + 1;
end

% testHigherOrderModes
ehb_high = ElegantHermiteBeam(w0, lambda, 2, 2);
field_high = ehb_high.opticalField(X, Y, 0.01);
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

% testZeroOrderElegant
ehb_00 = ElegantHermiteBeam(w0, lambda, 0, 0);
field_00 = ehb_00.opticalField(X, Y, 0);
if (all(all(isfinite(field_00))))
    fprintf('  PASS: zero order elegant valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: zero order elegant\n');
    failed = failed + 1;
end

% testGouyPhaseIncluded
params_z = ehb.getParameters(0.05);
if (params_z.GouyPhase > 0)
    fprintf('  PASS: Gouy phase included\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Gouy phase\n');
    failed = failed + 1;
end

% testDifferentNM
ehb_n2 = ElegantHermiteBeam(w0, lambda, 2, 0);
ehb_m2 = ElegantHermiteBeam(w0, lambda, 0, 2);
f_n2 = ehb_n2.opticalField(X, Y, 0);
f_m2 = ehb_m2.opticalField(X, Y, 0);
if (all(all(isfinite(f_n2))) && all(all(isfinite(f_m2))))
    fprintf('  PASS: different n m orders\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different n m\n');
    failed = failed + 1;
end

% testElegantAlphaFormula (alpha via ElegantHermiteParameters)
ehp_alpha = ElegantHermiteParameters(0.05, w0, lambda, 1, 1);
alpha_calc = ehp_alpha.alpha;
if (isfinite(alpha_calc))
    fprintf('  PASS: elegant alpha formula\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant alpha\n');
    failed = failed + 1;
end

% testWaistFromParameters
params_w = ehb.getParameters(0.1);
if (params_w.Waist > 0)
    fprintf('  PASS: waist from parameters\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist from parameters\n');
    failed = failed + 1;
end

% testElegantNegativeZ
field_neg = ehb.opticalField(X, Y, -0.05);
if (all(all(isfinite(field_neg))))
    fprintf('  PASS: elegant negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant negative z\n');
    failed = failed + 1;
end

% testElegantDifferentWavelengths
ehb_wl1 = ElegantHermiteBeam(w0, 532e-9, 1, 1);
ehb_wl2 = ElegantHermiteBeam(w0, 1064e-9, 1, 1);
f_wl1 = ehb_wl1.opticalField(X, Y, 0);
f_wl2 = ehb_wl2.opticalField(X, Y, 0);
if (all(all(isfinite(f_wl1))) && all(all(isfinite(f_wl2))))
    fprintf('  PASS: elegant different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant wavelengths\n');
    failed = failed + 1;
end

% testElegantHigherNM
ehb_hnm = ElegantHermiteBeam(w0, lambda, 3, 3);
f_hnm = ehb_hnm.opticalField(X, Y, 0);
if (all(all(isfinite(f_hnm))))
    fprintf('  PASS: elegant higher n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant higher n m\n');
    failed = failed + 1;
end

% testElegantFieldAmplitude
field_amp = ehb_00.opticalField(X, Y, 0);
max_amp_e = max(max(abs(field_amp)));
if (max_amp_e > 0)
    fprintf('  PASS: elegant field amplitude\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant amplitude\n');
    failed = failed + 1;
end

% testElegantAsymmetricNM
ehb_asym = ElegantHermiteBeam(w0, lambda, 2, 1);
f_asym = ehb_asym.opticalField(X, Y, 0);
if (all(all(isfinite(f_asym))))
    fprintf('  PASS: elegant asymmetric n m\n');
    passed = passed + 1;
else
    fprintf('  FAIL: elegant asymmetric\n');
    failed = failed + 1;
end

% testPhaseIncluded
field_phase = ehb.opticalField(X, Y, 0.05);
if (~isreal(field_phase) || any(any(abs(imag(field_phase)) > 0)))
    fprintf('  PASS: phase included in field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase included in field\n');
    failed = failed + 1;
end

% testFieldAtDifferentZ
f_z1 = ehb.opticalField(X, Y, 0.01);
f_z2 = ehb.opticalField(X, Y, 0.1);
if (all(all(isfinite(f_z1))) && all(all(isfinite(f_z2))))
    fprintf('  PASS: field at different z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field at different z\n');
    failed = failed + 1;
end

% testModalGouyRelativePhase (EHG10 vs EHG00)
z_rel = pi * w0^2 / lambda;
psi_rel = atan2(z_rel, pi * w0^2 / lambda);
x_rel = 0.7 * w0;
y_rel = 0;
ehb00_rel = ElegantHermiteBeam(w0, lambda, 0, 0);
ehb10_rel = ElegantHermiteBeam(w0, lambda, 1, 0);
f00_rel = ehb00_rel.opticalField(x_rel, y_rel, z_rel);
f10_rel = ehb10_rel.opticalField(x_rel, y_rel, z_rel);
k_rel = 2 * pi / lambda;
q_rel = z_rel + 1i * (pi * w0^2 / lambda);
alpha_rel = 1i * k_rel / (2 * q_rel);
H1_rel = PolynomialUtils.hermitePoly(1, sqrt(alpha_rel) * x_rel);
ratio_rel = f10_rel / (H1_rel * f00_rel);
phase_err_rel = angle(ratio_rel * exp(1i * psi_rel));
if (abs(phase_err_rel) < 1e-8)
    fprintf('  PASS: modal Gouy relative phase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: modal Gouy relative phase\n');
    failed = failed + 1;
end

% testBeamName
if (strcmp(ehb.beamName(), 'elegant_hermite_1_1'))
    fprintf('  PASS: beamName\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beamName\n');
    failed = failed + 1;
end

fprintf('\n=== ElegantHermiteBeam: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
