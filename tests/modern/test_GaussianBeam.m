% Compatible with GNU Octave and MATLAB
% Tests for GaussianBeam (Phase 3 API: GaussianBeam(w0, lambda))

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== GaussianBeam Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

beam = GaussianBeam(w0, lambda);

% testFieldGeneration
field = beam.opticalField(X, Y, 0.01);
if (size(field) == [64, 64])
    fprintf('  PASS: field generation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field generation\n');
    failed = failed + 1;
end

% testAmplitudeAtCenter
field0 = beam.opticalField(zeros(64,64), zeros(64,64), 0);
if (abs(abs(field0(33,33)) - 1) < 1e-10)
    fprintf('  PASS: amplitude at center\n');
    passed = passed + 1;
else
    fprintf('  FAIL: amplitude at center\n');
    failed = failed + 1;
end

% testDecaysWithRadius
X_far = 3*w0*ones(64,64);
field_far = beam.opticalField(X_far, zeros(64,64), 0);
if (abs(field_far(1,1)) < 0.01)
    fprintf('  PASS: decays with radius\n');
    passed = passed + 1;
else
    fprintf('  FAIL: decays with radius\n');
    failed = failed + 1;
end

% testAtZ0Finite
field_z0 = beam.opticalField(X, Y, 0);
if (all(all(isfinite(field_z0))))
    fprintf('  PASS: at z=0 finite\n');
    passed = passed + 1;
else
    fprintf('  FAIL: at z=0 finite\n');
    failed = failed + 1;
end

% testInitialWaistStored
if (beam.InitialWaist == w0)
    fprintf('  PASS: InitialWaist stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: InitialWaist stored\n');
    failed = failed + 1;
end

% testComplexField
field_c = beam.opticalField(X, Y, 0.05);
if (~isreal(field_c))
    fprintf('  PASS: complex field\n');
    passed = passed + 1;
else
    fprintf('  FAIL: complex field\n');
    failed = failed + 1;
end

% testPhaseAtWaist
field_p = beam.opticalField(zeros(64,64), zeros(64,64), 0);
phase_center = angle(field_p(33,33));
if (abs(phase_center) < pi/2)
    fprintf('  PASS: phase at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: phase at waist\n');
    failed = failed + 1;
end

% testGetParameters
params = beam.getParameters(0.1);
if (abs(params.zCoordinate - 0.1) < 1e-15 && params.InitialWaist == w0)
    fprintf('  PASS: getParameters(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getParameters(z)\n');
    failed = failed + 1;
end

% testFieldSymmetry
field_sym = beam.opticalField(X, Y, 0);
field_center = abs(field_sym(33,33));
field_edge   = abs(field_sym(1,1));
if (field_edge < field_center)
    fprintf('  PASS: field symmetry\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field symmetry\n');
    failed = failed + 1;
end

% testWaistPropagates
params_w1 = beam.getParameters(0);
params_w2 = beam.getParameters(0.1);
if (params_w2.Waist > params_w1.Waist)
    fprintf('  PASS: waist propagates\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waist propagates\n');
    failed = failed + 1;
end

% testDifferentWavelength
beam_l2 = GaussianBeam(w0, 1064e-9);
field_l2 = beam_l2.opticalField(X, Y, 0);
if (all(all(isfinite(field_l2))))
    fprintf('  PASS: different wavelength\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different wavelength\n');
    failed = failed + 1;
end

% testDifferentWaist — smaller w0 → smaller Rayleigh distance
beam_w02 = GaussianBeam(50e-6, lambda);
params_ref = beam.getParameters(0);
params_w02 = beam_w02.getParameters(0);
if (params_w02.RayleighDistance < params_ref.RayleighDistance)
    fprintf('  PASS: different waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different waist\n');
    failed = failed + 1;
end

% testFieldWithPositiveZ
field_pos = beam.opticalField(X, Y, 0.05);
if (all(all(isfinite(field_pos))))
    fprintf('  PASS: field with positive z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: positive z\n');
    failed = failed + 1;
end

% testFieldWithNegativeZ
field_neg = beam.opticalField(X, Y, -0.05);
if (all(all(isfinite(field_neg))))
    fprintf('  PASS: field with negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative z\n');
    failed = failed + 1;
end

% testConstructorEmpty
try
    beam_empty = GaussianBeam();
    fprintf('  PASS: constructor empty\n');
    passed = passed + 1;
catch
    fprintf('  FAIL: constructor empty\n');
    failed = failed + 1;
end

% testBeamDivergence
X_div = 5*w0*ones(64,64);
field_div = beam.opticalField(X_div, zeros(64,64), 0);
if (abs(field_div(1,1)) < 0.1)
    fprintf('  PASS: beam divergence\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beam divergence\n');
    failed = failed + 1;
end

% testFieldNormalized
field_norm = beam.opticalField(zeros(64,64), zeros(64,64), 0);
max_amp = max(max(abs(field_norm)));
if (abs(max_amp - 1) < 1e-10)
    fprintf('  PASS: field normalized\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field normalized\n');
    failed = failed + 1;
end

% testBeamName
if (strcmp(beam.beamName(), 'gaussian'))
    fprintf('  PASS: beamName\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beamName\n');
    failed = failed + 1;
end

fprintf('\n=== GaussianBeam: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
