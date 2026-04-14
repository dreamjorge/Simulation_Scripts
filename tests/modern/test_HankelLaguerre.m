% Compatible with GNU Octave and MATLAB
% Tests for HankelLaguerre (Phase 3 API: HankelLaguerre(w0, lambda, l, p, hankelType))

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

fprintf('=== HankelLaguerre Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 1e-3;
lambda = 632e-9;
grid = GridUtils(64, 64, 1e-3, 1e-3);
[X, Y] = grid.create2DGrid();

% testConstructorHankelType1
HL1 = HankelLaguerre(w0, lambda, 1, 0, 1);
field1 = HL1.opticalField(X, Y, 0);
if (~isempty(field1) && size(field1, 1) == 64 && size(field1, 2) == 64)
    fprintf('  PASS: constructor hankelType=1\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor hankelType=1\n');
    failed = failed + 1;
end

% testConstructorHankelType2
HL2 = HankelLaguerre(w0, lambda, 1, 0, 2);
field2 = HL2.opticalField(X, Y, 0);
if (~isempty(field2) && size(field2) == [64, 64])
    fprintf('  PASS: constructor hankelType=2\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor hankelType=2\n');
    failed = failed + 1;
end

% testConstructorDefaultType
HL_def = HankelLaguerre(w0, lambda, 1, 0);
field_def = HL_def.opticalField(X, Y, 0);
if (~isempty(field_def))
    fprintf('  PASS: constructor default type\n');
    passed = passed + 1;
else
    fprintf('  FAIL: default type failed\n');
    failed = failed + 1;
end

% testHankelTypeStored
if (HL1.HankelType == 1 && HL2.HankelType == 2)
    fprintf('  PASS: HankelType stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelType stored\n');
    failed = failed + 1;
end

% testModeIndicesStored
if (HL1.l == 1 && HL1.p == 0)
    fprintf('  PASS: mode indices stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: mode indices stored\n');
    failed = failed + 1;
end

% testInitialWaistStored
if (HL1.InitialWaist == w0)
    fprintf('  PASS: InitialWaist stored\n');
    passed = passed + 1;
else
    fprintf('  FAIL: InitialWaist stored\n');
    failed = failed + 1;
end

% testGetParameters
params = HL1.getParameters(0.1);
if (abs(params.zCoordinate - 0.1) < 1e-15 && params.InitialWaist == w0)
    fprintf('  PASS: getParameters(z)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: getParameters(z)\n');
    failed = failed + 1;
end

% testHankelType1PlusType2 (H1 + H2 = 2*LG)
% H^(1) = LG + i*XLG,  H^(2) = LG - i*XLG  =>  H1 + H2 = 2*LG
lb = LaguerreBeam(w0, lambda, 1, 0);
LG_field = lb.opticalField(X, Y, 0);
sum_field = field1 + field2;
if (max(max(abs(sum_field - 2*LG_field))) < 1e-10)
    fprintf('  PASS: H1 + H2 = 2*LG\n');
    passed = passed + 1;
else
    fprintf('  FAIL: H1 + H2 != 2*LG\n');
    failed = failed + 1;
end

% testFieldFiniteAtWaist
f_w1 = HL1.opticalField(X, Y, 0);
f_w2 = HL2.opticalField(X, Y, 0);
if (all(all(isfinite(f_w1))) && all(all(isfinite(f_w2))))
    fprintf('  PASS: finite at waist\n');
    passed = passed + 1;
else
    fprintf('  FAIL: finite at waist\n');
    failed = failed + 1;
end

% testFieldFiniteAtPropagation
f_p1 = HL1.opticalField(X, Y, 0.1);
f_p2 = HL2.opticalField(X, Y, 0.1);
if (all(all(isfinite(f_p1))) && all(all(isfinite(f_p2))))
    fprintf('  PASS: finite at propagation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: finite at propagation\n');
    failed = failed + 1;
end

% testHigherModes
HL_high = HankelLaguerre(w0, lambda, 2, 1, 1);
f_high = HL_high.opticalField(X, Y, 0);
if (all(all(isfinite(f_high))))
    fprintf('  PASS: higher modes valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: higher modes\n');
    failed = failed + 1;
end

% testNegativeL
HL_neg = HankelLaguerre(w0, lambda, -1, 0, 1);
f_neg = HL_neg.opticalField(X, Y, 0);
if (all(all(isfinite(f_neg))))
    fprintf('  PASS: negative l valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: negative l\n');
    failed = failed + 1;
end

% testL0P0
HL_00 = HankelLaguerre(w0, lambda, 0, 0, 1);
f_00 = HL_00.opticalField(X, Y, 0);
if (all(all(isfinite(f_00))))
    fprintf('  PASS: l=0 p=0 valid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: l=0 p=0\n');
    failed = failed + 1;
end

% testFieldAmplitude
max_amp = max(max(abs(field1)));
if (max_amp > 0)
    fprintf('  PASS: field amplitude positive\n');
    passed = passed + 1;
else
    fprintf('  FAIL: field amplitude\n');
    failed = failed + 1;
end

% testDifferentWavelengths
HL_lam1 = HankelLaguerre(w0, 532e-9, 1, 0, 1);
HL_lam2 = HankelLaguerre(w0, 1064e-9, 1, 0, 1);
f_lam1 = HL_lam1.opticalField(X, Y, 0);
f_lam2 = HL_lam2.opticalField(X, Y, 0);
if (all(all(isfinite(f_lam1))) && all(all(isfinite(f_lam2))))
    fprintf('  PASS: different wavelengths\n');
    passed = passed + 1;
else
    fprintf('  FAIL: different wavelengths\n');
    failed = failed + 1;
end

% testBeamName
if (strcmp(HL1.beamName(), 'hankel1_laguerre_1_0'))
    fprintf('  PASS: beamName\n');
    passed = passed + 1;
else
    fprintf('  FAIL: beamName\n');
    failed = failed + 1;
end

fprintf('\n=== HankelLaguerre: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
