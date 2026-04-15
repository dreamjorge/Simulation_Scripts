% Compatible with GNU Octave and MATLAB
% Tests for GridUtils

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

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
if (isequal(size(Kx), [Nx, Nx]))
    fprintf('  PASS: createFreqGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: createFreqGrid\n');
    failed = failed + 1;
end

% testMeshgrid2D
[Xm, Ym] = GridUtils.meshgrid2D(32, 1e-3);
if (isequal(size(Xm), [32, 32]))
    fprintf('  PASS: meshgrid2D\n');
    passed = passed + 1;
else
    fprintf('  FAIL: meshgrid2D\n');
    failed = failed + 1;
end

% testFreqGrid
[Kxf, Kyf] = GridUtils.freqGrid(32, 1e-3);
if (isequal(size(Kxf), [32, 32]))
    fprintf('  PASS: freqGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freqGrid\n');
    failed = failed + 1;
end

% testPolarGrid
[r, theta] = GridUtils.polarGrid(32, 1e-3);
if (isequal(size(r), [32, 32]))
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

% testNxNyProperties
if (grid.Nx == Nx && grid.Ny == Ny)
    fprintf('  PASS: Nx Ny properties\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Nx Ny\n');
    failed = failed + 1;
end

% testDxDyProperties
if (grid.Dx == Dx && grid.Dy == Dy)
    fprintf('  PASS: Dx Dy properties\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Dx Dy\n');
    failed = failed + 1;
end

% testDxDyCalculated
if (grid.dx > 0 && grid.dy > 0)
    fprintf('  PASS: dx dy calculated\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dx dy\n');
    failed = failed + 1;
end

% testNonSquareGrid
grid_ns = GridUtils(128, 64, 2e-3, 1e-3);
[Xns, Yns] = grid_ns.create2DGrid();
if (size(Xns,1) == 64 && size(Xns,2) == 128)
    fprintf('  PASS: non-square grid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: non-square\n');
    failed = failed + 1;
end

% testCreate2DGridCenterAtZero
grid_c = GridUtils(64, 64, 1e-3, 1e-3);
[Xc, Yc] = grid_c.create2DGrid();
if (abs(Xc(33,33)) < 1e-10 && abs(Yc(33,33)) < 1e-10)
    fprintf('  PASS: center at zero\n');
    passed = passed + 1;
else
    fprintf('  FAIL: center\n');
    failed = failed + 1;
end

% testFreqGridAtOrigin
grid_f = GridUtils(64, 64, 1e-3, 1e-3);
[Kxf, Kyf] = grid_f.createFreqGrid();
if (abs(Kxf(33,33)) < 1 && abs(Kyf(33,33)) < 1)
    fprintf('  PASS: freq grid origin near zero\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freq origin\n');
    failed = failed + 1;
end

% testPolarGridCenter
[rpc, tpc] = GridUtils.polarGrid(64, 1e-3);
if (rpc(33,33) == 0 && tpc(33,33) == 0)
    fprintf('  PASS: polar grid center at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: polar center\n');
    failed = failed + 1;
end

% test3DGridDz
grid_3d = GridUtils(16, 16, 1e-3, 1e-3, 8, 2e-3);
if (grid_3d.dz > 0)
    fprintf('  PASS: 3D grid dz\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 3D dz\n');
    failed = failed + 1;
end

% testCreate2DGridValues
grid_v = GridUtils(8, 8, 1e-3, 1e-3);
[Xv, Yv] = grid_v.create2DGrid();
if (Xv(1,1) < 0 && Xv(1,8) > 0 && Yv(1,1) < 0 && Yv(8,1) > 0)
    fprintf('  PASS: create2DGrid values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 2D values\n');
    failed = failed + 1;
end

% testCreateFreqGridValues
grid_fv = GridUtils(8, 8, 1e-3, 1e-3);
[Kxfv, Kyfv] = grid_fv.createFreqGrid();
if (Kxfv(5,5) == 0 && Kyfv(5,5) == 0)
    fprintf('  PASS: createFreqGrid values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freq values\n');
    failed = failed + 1;
end

% testMeshgrid2DStaticValues
[Xms, Yms] = GridUtils.meshgrid2D(4, 1e-3);
if (isequal(size(Xms), [4, 4]) && isequal(size(Yms), [4, 4]))
    fprintf('  PASS: meshgrid2D static\n');
    passed = passed + 1;
else
    fprintf('  FAIL: meshgrid2D static\n');
    failed = failed + 1;
end

% testFreqGridStaticValues
[Kxs, Kys] = GridUtils.freqGrid(4, 1e-3);
if (Kxs(3,3) == 0 && Kys(3,3) == 0)
    fprintf('  PASS: freqGrid static values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freqGrid static\n');
    failed = failed + 1;
end

% testPolarGridValues
[rpv, tpv] = GridUtils.polarGrid(8, 1e-3);
if (rpv(5,5) == 0 && abs(tpv(5,5)) < 1e-10)
    fprintf('  PASS: polarGrid values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: polar values\n');
    failed = failed + 1;
end

% testGridWithoutNz
grid_no3d = GridUtils(16, 16, 1e-3, 1e-3);
if (~isfield(grid_no3d, 'Nz') || grid_no3d.Nz == 0)
    fprintf('  PASS: grid without Nz\n');
    passed = passed + 1;
else
    fprintf('  FAIL: no Nz\n');
    failed = failed + 1;
end

% testDifferentDxDyManual
dx_manual = 1e-3 / 16;
dy_manual = 1e-3 / 32;
if (abs(dx_manual - dy_manual) > 1e-10)
    fprintf('  PASS: different dx dy\n');
    passed = passed + 1;
else
    fprintf('  FAIL: dx dy\n');
    failed = failed + 1;
end

% test3DGridSize
grid_3ds = GridUtils(8, 8, 1e-3, 1e-3, 4, 1e-3);
[X3s, Y3s, Z3s] = grid_3ds.create3DGrid();
if (size(X3s,1) == 8 && size(X3s,2) == 8 && size(Z3s,3) == 4)
    fprintf('  PASS: 3D grid size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 3D size\n');
    failed = failed + 1;
end

% test3DGridZValues
grid_zv = GridUtils(8, 8, 1e-3, 1e-3, 4, 4e-3);
[~, ~, Zzv] = grid_zv.create3DGrid();
if (Zzv(1,1,1) == 0 && Zzv(1,1,4) > Zzv(1,1,1))
    fprintf('  PASS: 3D grid z values\n');
    passed = passed + 1;
else
    fprintf('  FAIL: z values\n');
    failed = failed + 1;
end

% testStaticMeshgrid2DDifferentN
[Xmd, Ymd] = GridUtils.meshgrid2D(16, 2e-3);
if (isequal(size(Xmd), [16, 16]))
    fprintf('  PASS: meshgrid2D different N\n');
    passed = passed + 1;
else
    fprintf('  FAIL: meshgrid2D N\n');
    failed = failed + 1;
end

% testStaticFreqGridDifferentN
[Kxd, Kyd] = GridUtils.freqGrid(16, 2e-3);
if (isequal(size(Kxd), [16, 16]))
    fprintf('  PASS: freqGrid different N\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freqGrid N\n');
    failed = failed + 1;
end

% testCreate2DGridSymmetry
grid_s = GridUtils(8, 8, 1e-3, 1e-3);
[Xs, Ys] = grid_s.create2DGrid();
if (min(min(Xs)) < 0 && max(max(Xs)) > 0 && min(min(Ys)) < 0 && max(max(Ys)) > 0)
    fprintf('  PASS: create2DGrid symmetry\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 2D grid symmetry\n');
    failed = failed + 1;
end

% testCreateFreqGridSymmetry
grid_fs = GridUtils(8, 8, 1e-3, 1e-3);
[Kxs_f, Kys_f] = grid_fs.createFreqGrid();
if (min(min(Kxs_f)) < 0 && max(max(Kxs_f)) > 0 && min(min(Kys_f)) < 0 && max(max(Kys_f)) > 0)
    fprintf('  PASS: createFreqGrid symmetry\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freq grid symmetry\n');
    failed = failed + 1;
end

% testPolarGridSymmetry
[rps, tps] = GridUtils.polarGrid(8, 1e-3);
if (max(max(rps)) > 0)
    fprintf('  PASS: polar grid symmetry\n');
    passed = passed + 1;
else
    fprintf('  FAIL: polar symmetry\n');
    failed = failed + 1;
end

% testGridPropertiesNx
grid_p = GridUtils(32, 64, 2e-3, 4e-3);
if (grid_p.Nx == 32 && grid_p.Ny == 64)
    fprintf('  PASS: grid properties Nx Ny\n');
    passed = passed + 1;
else
    fprintf('  FAIL: grid properties\n');
    failed = failed + 1;
end

% testGridPropertiesDxDy
if (grid_p.Dx == 2e-3 && grid_p.Dy == 4e-3)
    fprintf('  PASS: grid properties Dx Dy\n');
    passed = passed + 1;
else
    fprintf('  FAIL: Dx Dy properties\n');
    failed = failed + 1;
end

% testGridCalculatedDxDy
if (abs(grid_p.dx - 2e-3/32) < 1e-10 && abs(grid_p.dy - 4e-3/64) < 1e-10)
    fprintf('  PASS: grid calculated dx dy\n');
    passed = passed + 1;
else
    fprintf('  FAIL: calculated dx dy\n');
    failed = failed + 1;
end

% test3DGridWithLargeNz
grid_3dl = GridUtils(8, 8, 1e-3, 1e-3, 100, 1e-2);
if (grid_3dl.Nz == 100)
    fprintf('  PASS: 3D grid large Nz\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 3D large Nz\n');
    failed = failed + 1;
end

% test3DGridDzCalculated
if (abs(grid_3dl.dz - 1e-2/100) < 1e-15)
    fprintf('  PASS: 3D grid dz calculated\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 3D dz calc\n');
    failed = failed + 1;
end

% testMeshgrid2DValuesRange
[Xmr, Ymr] = GridUtils.meshgrid2D(4, 1e-3);
if (min(min(Xmr)) < 0 && max(max(Xmr)) > 0 && min(min(Ymr)) < 0 && max(max(Ymr)) > 0)
    fprintf('  PASS: meshgrid2D values range\n');
    passed = passed + 1;
else
    fprintf('  FAIL: meshgrid2D range\n');
    failed = failed + 1;
end

% testFreqGridValuesRange
[Kxr, Kyr] = GridUtils.freqGrid(4, 1e-3);
if (min(min(Kxr)) < 0 && max(max(Kxr)) > 0 && min(min(Kyr)) < 0 && max(max(Kyr)) > 0)
    fprintf('  PASS: freqGrid values range\n');
    passed = passed + 1;
else
    fprintf('  FAIL: freqGrid range\n');
    failed = failed + 1;
end

% testCreate2DGridSinglePoint
grid_sp = GridUtils(2, 2, 1e-6, 1e-6);
[Xsp, Ysp] = grid_sp.create2DGrid();
if (numel(Xsp) == 4)
    fprintf('  PASS: create2DGrid single point\n');
    passed = passed + 1;
else
    fprintf('  FAIL: 2D single point\n');
    failed = failed + 1;
end

% testGridWithVerySmallDx
grid_vs = GridUtils(64, 64, 1e-9, 1e-9);
if (grid_vs.dx > 0)
    fprintf('  PASS: grid very small dx\n');
    passed = passed + 1;
else
    fprintf('  FAIL: very small dx\n');
    failed = failed + 1;
end

% testAsymmetricMeshgrid2D
[Xa, Ya] = GridUtils.meshgrid2D(32, 64, 1e-3, 2e-3);
if (size(Xa,1) == 64 && size(Xa,2) == 32 && size(Ya,1) == 64 && size(Ya,2) == 32)
    fprintf('  PASS: asymmetric meshgrid2D\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric meshgrid2D\n');
    failed = failed + 1;
end

% testAsymmetricFreqGrid
[Kxa, Kya] = GridUtils.freqGrid(32, 64, 1e-3, 2e-3);
if (size(Kxa,1) == 64 && size(Kxa,2) == 32 && size(Kya,1) == 64 && size(Kya,2) == 32)
    fprintf('  PASS: asymmetric freqGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric freqGrid\n');
    failed = failed + 1;
end

% testAsymmetricPolarGrid
[rpa, tpa] = GridUtils.polarGrid(32, 64, 1e-3, 2e-3);
if (size(rpa,1) == 64 && size(rpa,2) == 32)
    fprintf('  PASS: asymmetric polarGrid\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric polarGrid\n');
    failed = failed + 1;
end

% testAsymmetricMeshgrid2D_3args
[Xa3, Ya3] = GridUtils.meshgrid2D(32, 64, 1e-3);
if (size(Xa3,1) == 64 && size(Xa3,2) == 32)
    fprintf('  PASS: asymmetric meshgrid2D 3 args\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric meshgrid2D 3 args\n');
    failed = failed + 1;
end

% testAsymmetricFreqGrid_3args
[Kxa3, Kya3] = GridUtils.freqGrid(32, 64, 1e-3);
if (size(Kxa3,1) == 64 && size(Kxa3,2) == 32)
    fprintf('  PASS: asymmetric freqGrid 3 args\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric freqGrid 3 args\n');
    failed = failed + 1;
end

% testAsymmetricPolarGrid_3args
[rpa3, tpa3] = GridUtils.polarGrid(32, 64, 1e-3);
if (size(rpa3,1) == 64 && size(rpa3,2) == 32)
    fprintf('  PASS: asymmetric polarGrid 3 args\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric polarGrid 3 args\n');
    failed = failed + 1;
end

% testAsymmetricGridValues
[Xav, Yav] = GridUtils.meshgrid2D(4, 8, 1e-3, 2e-3);
x_range_ok = (min(min(Xav)) < 0) && (max(max(Xav)) > 0);
y_range_ok = (min(min(Yav)) < 0) && (max(max(Yav)) > 0);
if x_range_ok && y_range_ok
    fprintf('  PASS: asymmetric grid values range\n');
    passed = passed + 1;
else
    fprintf('  FAIL: asymmetric grid values range\n');
    failed = failed + 1;
end

fprintf('\n=== GridUtils: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
