#!/usr/bin/env octave
% Tests for HankelLaguerre

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== HankelLaguerre Tests ===\n\n');
passed = 0;
failed = 0;

N = 64; D = 1e-3;
grid = GridUtils(N, N, D, D);
[X, Y] = grid.create2DGrid();
[r, theta] = cart2pol(X, Y);

params = LaguerreParameters(0, 1e-3, 632e-9, 1, 0);

% testConstructorHankelType1
try
    HL = HankelLaguerre(r, theta, params, 1);
    if ~isempty(HL.OpticalFieldLaguerre) && size(HL.OpticalFieldLaguerre) == size(r)
        fprintf('  PASS: constructor hankelType=1\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: OpticalFieldLaguerre empty or wrong size\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testConstructorHankelType2
try
    HL = HankelLaguerre(r, theta, params, 2);
    if ~isempty(HL.OpticalFieldLaguerre)
        fprintf('  PASS: constructor hankelType=2\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: OpticalFieldLaguerre empty\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testConstructorDefaultType
try
    HL = HankelLaguerre(r, theta, params);
    if ~isempty(HL.OpticalFieldLaguerre)
        fprintf('  PASS: constructor default type\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: default type failed\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

% testHankelType1PlusType2
try
    HL1 = HankelLaguerre(r, theta, params, 1);
    HL2 = HankelLaguerre(r, theta, params, 2);
    % H1 + H2 should equal 2*LB (imag parts cancel)
    diff = HL1.OpticalFieldLaguerre + HL2.OpticalFieldLaguerre;
    LB = LaguerreBeam(r, theta, params);
    if max(max(abs(diff - 2*LB.OpticalField))) < 1e-10
        fprintf('  PASS: H1 + H2 = 2*LB\n');
        passed = passed + 1;
    else
        fprintf('  FAIL: H1 + H2 != 2*LB\n');
        failed = failed + 1;
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    failed = failed + 1;
end

fprintf('\n=== HankelLaguerre: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end