#!/usr/bin/env octave
% Tests for GridUtils

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== GridUtils Tests ===\n\n');
passed = 0;
failed = 0;

Nx = 64; Ny = 64; Dx = 1e-3; Dy = 1e-3;
grid = GridUtils(Nx, Ny, Dx, Dy);

% testCreate2DGrid
[X, Y] = grid.create2DGrid();
if (size(X,1) == Ny && size(X,2) == Nx)
    fprintf('  PASS: create2DGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: create2DGrid\n');
    failed = failed + 1;
end

% testCreateFreqGrid
[Kx, Ky] = grid.createFreqGrid();
if (size(Kx) == [Nx, Nx])
    fprintf('  PASS: createFreqGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: createFreqGrid\n');
    failed = failed + 1;
end

% testMeshgrid2D
[Xm, Ym] = GridUtils.meshgrid2D(32, 1e-3);
if (size(Xm) == [32, 32])
    fprintf('  PASS: meshgrid2D\n');
    passed = passed + 1;
else
    fprintf('  FAIL: meshgrid2D\n');
    failed = failed + 1;
end

% testFreqGrid
[Kxf, Kyf] = GridUtils.freqGrid(32, 1e-3);
if (size(Kxf) == [32, 32])
    fprintf('  PASS: freqGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freqGrid\n');
    failed = failed + 1;
end

% testPolarGrid
[r, theta] = GridUtils.polarGrid(32, 1e-3);
if (size(r) == [32, 32])
    fprintf('  PASS: polarGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: polarGrid\n');
    failed = failed + 1;
end

% testCreate3DGrid
grid3d = GridUtils(16, 16, 1e-3, 1e-3, 8, 1e-3);
[X3, Y3, Z3] = grid3d.create3DGrid();
if (ndims(X3) == 3)
    fprintf('  PASS: create3DGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: create3DGrid\n');
    failed = failed + 1;
end

fprintf('\n=== GridUtils: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
