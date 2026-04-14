% Compatible with GNU Octave and MATLAB
% Legacy compatibility tests for Hankel Hermite/Laguerre wrappers.

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== Hankel Compatibility Tests ===\n\n');
passed = 0;
failed = 0;

w0 = 100e-6;
lambda = 632.8e-9;
z = 0.01;

% testHermiteLegacySolutions
[h0, nh0] = HermiteParameters.getHermiteSolutions(2, linspace(-1, 1, 11));
if (~isempty(h0) && ~isempty(nh0) && numel(h0) == 11 && numel(nh0) == 11)
    fprintf('  PASS: HermiteParameters.getHermiteSolutions\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HermiteParameters.getHermiteSolutions\n');
    failed = failed + 1;
end

% testHankelHermiteConstructor
hp = HermiteParameters(z, w0, lambda, 1, 1);
x = linspace(-3e-4, 3e-4, 51);
y = 1e-5;
hh11 = HankelHermite(x, y, hp, 11);
if (~isempty(hh11.OpticalField) && numel(hh11.OpticalField) == numel(x))
    fprintf('  PASS: HankelHermite constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelHermite constructor\n');
    failed = failed + 1;
end

% testHankeleHermiteAlias
hhe11 = HankeleHermite(x, y, hp, 11);
if (max(abs(hhe11.OpticalField - hh11.OpticalField)) < 1e-12)
    fprintf('  PASS: HankeleHermite alias\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankeleHermite alias\n');
    failed = failed + 1;
end

% testHankelLaguerreLegacyConstructor
lp = LaguerreParameters(z, w0, lambda, 1, 0);
r = linspace(0, 3e-4, 51);
th = pi / 4;
hl1 = HankelLaguerre(r, th, lp, 1);
if (~isempty(hl1.OpticalFieldLaguerre) && numel(hl1.OpticalFieldLaguerre) == numel(r))
    fprintf('  PASS: HankelLaguerre legacy constructor\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelLaguerre legacy constructor\n');
    failed = failed + 1;
end

% testHankeleLaguerreAlias
hle1 = HankeleLaguerre(r, th, lp, 1);
if (max(abs(hle1.OpticalFieldLaguerre - hl1.OpticalFieldLaguerre)) < 1e-12)
    fprintf('  PASS: HankeleLaguerre alias\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankeleLaguerre alias\n');
    failed = failed + 1;
end

fprintf('\n=== Hankel Compatibility: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
