% Compatible with GNU Octave and MATLAB
% Tests for PhysicalConstants

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'ParaxialBeams'));

fprintf('=== PhysicalConstants Tests ===\n\n');
passed = 0;
failed = 0;

lambda = 632.8e-9;
w0 = 100e-6;
zr = pi * w0^2 / lambda;

% testWaveNumber
k = PhysicalConstants.waveNumber(lambda);
expected = 2*pi / lambda;
if (abs(k - expected) < 1e-5)
    fprintf('  PASS: waveNumber\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber\n');
    failed = failed + 1;
end

% testRayleighDistance
zr_test = PhysicalConstants.rayleighDistance(w0, lambda);
expected_zr = pi * w0^2 / lambda;
if (abs(zr_test - expected_zr) < 1e-10)
    fprintf('  PASS: rayleighDistance\n');
    passed = passed + 1;
else
    fprintf('  FAIL: rayleighDistance\n');
    failed = failed + 1;
end

% testWaistAtZ
w = PhysicalConstants.waistAtZ(w0, 0.05, lambda);
expected_w = w0 * sqrt(1 + (0.05/zr)^2);
if (abs(w - expected_w) / expected_w < 1e-5)
    fprintf('  PASS: waistAtZ\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ\n');
    failed = failed + 1;
end

% testRadiusOfCurvature
R = PhysicalConstants.radiusOfCurvature(0.05, zr);
if (R > 0)
    fprintf('  PASS: radiusOfCurvature\n');
    passed = passed + 1;
else
    fprintf('  FAIL: radiusOfCurvature\n');
    failed = failed + 1;
end

% testGouyPhase
gouy = PhysicalConstants.gouyPhase(0.05, zr);
if (gouy > 0)
    fprintf('  PASS: gouyPhase\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase\n');
    failed = failed + 1;
end

% test speed_of_light
c = PhysicalConstants.speed_of_light;
if (c > 1e8 && c < 3e8)
    fprintf('  PASS: speed_of_light\n');
    passed = passed + 1;
else
    fprintf('  FAIL: speed_of_light\n');
    failed = failed + 1;
end

% test planck
h = PhysicalConstants.planck;
if (h > 1e-34 && h < 1e-33)
    fprintf('  PASS: planck\n');
    passed = passed + 1;
else
    fprintf('  FAIL: planck\n');
    failed = failed + 1;
end

% test planck_reduced
hbar = PhysicalConstants.planck_reduced;
if (hbar > 1e-34 && hbar < 1e-33)
    fprintf('  PASS: planck_reduced\n');
    passed = passed + 1;
else
    fprintf('  FAIL: planck_reduced\n');
    failed = failed + 1;
end

% test vacuum_permittivity
eps0 = PhysicalConstants.vacuum_permittivity;
if (eps0 > 8e-12 && eps0 < 9e-12)
    fprintf('  PASS: vacuum_permittivity\n');
    passed = passed + 1;
else
    fprintf('  FAIL: vacuum_permittivity\n');
    failed = failed + 1;
end

% test vacuum_permeability
mu0 = PhysicalConstants.vacuum_permeability;
if (mu0 > 1e-6 && mu0 < 2e-6)
    fprintf('  PASS: vacuum_permeability\n');
    passed = passed + 1;
else
    fprintf('  FAIL: vacuum_permeability\n');
    failed = failed + 1;
end

% test impedance_vacuum
eta0 = PhysicalConstants.impedance_vacuum;
if (eta0 > 376 && eta0 < 377)
    fprintf('  PASS: impedance_vacuum\n');
    passed = passed + 1;
else
    fprintf('  FAIL: impedance_vacuum\n');
    failed = failed + 1;
end

% testWaistAtZWithoutZr
w_no_zr = PhysicalConstants.waistAtZ(w0, 0.05, lambda);
if (w_no_zr > w0)
    fprintf('  PASS: waistAtZ without zr\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ without zr\n');
    failed = failed + 1;
end

% testWaistAtZAtOrigin
w0_test = PhysicalConstants.waistAtZ(w0, 0, lambda);
if (abs(w0_test - w0) / w0 < 1e-10)
    fprintf('  PASS: waistAtZ at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ at origin\n');
    failed = failed + 1;
end

% testRadiusOfCurvatureAtOrigin
R0 = PhysicalConstants.radiusOfCurvature(0, zr);
if (isinf(R0))
    fprintf('  PASS: radiusOfCurvature at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: radiusOfCurvature at origin\n');
    failed = failed + 1;
end

% testGouyPhaseAtOrigin
gouy0 = PhysicalConstants.gouyPhase(0, zr);
if (abs(gouy0) < 1e-10)
    fprintf('  PASS: gouyPhase at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase at origin\n');
    failed = failed + 1;
end

% testRayleighDistanceVector
w0_vec = [50e-6, 100e-6, 200e-6];
zr_vec = PhysicalConstants.rayleighDistance(w0_vec, lambda);
if (numel(zr_vec) == 3 && all(zr_vec > 0))
    fprintf('  PASS: rayleighDistance vector\n');
    passed = passed + 1;
else
    fprintf('  FAIL: rayleighDistance vector\n');
    failed = failed + 1;
end

% testWaveNumberVector
lambda_vec = [532e-9, 632.8e-9, 1064e-9];
k_vec = PhysicalConstants.waveNumber(lambda_vec);
if (numel(k_vec) == 3 && all(k_vec > 0))
    fprintf('  PASS: waveNumber vector\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber vector\n');
    failed = failed + 1;
end

% testGouyPhaseNegativeZ
gouy_neg = PhysicalConstants.gouyPhase(-0.05, zr);
if (gouy_neg < 0)
    fprintf('  PASS: gouyPhase negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase negative z\n');
    failed = failed + 1;
end

% testRadiusOfCurvatureNegativeZ
R_neg = PhysicalConstants.radiusOfCurvature(-0.05, zr);
if (R_neg < 0)
    fprintf('  PASS: radiusOfCurvature negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: radiusOfCurvature negative z\n');
    failed = failed + 1;
end

% testWaistAtZVector
z_vec = [0, 0.05, 0.1, 0.2];
w_vec = PhysicalConstants.waistAtZ(w0, z_vec, lambda);
if (numel(w_vec) == 4 && all(w_vec >= w0))
    fprintf('  PASS: waistAtZ vector z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ vector z\n');
    failed = failed + 1;
end

% testWaistAtZVeryLargeZ
w_vl = PhysicalConstants.waistAtZ(w0, 10, lambda);
if (w_vl > w0 * sqrt(1 + (10/zr)^2) * 0.99)
    fprintf('  PASS: waistAtZ very large z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ very large z\n');
    failed = failed + 1;
end

% testGouyPhaseLargeZ
gouy_l = PhysicalConstants.gouyPhase(10*zr, zr);
if (gouy_l > 0 && gouy_l <= pi/2)
    fprintf('  PASS: gouyPhase large z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase large z\n');
    failed = failed + 1;
end

% testGouyPhaseVerySmallZ
gouy_vs = PhysicalConstants.gouyPhase(1e-10, zr);
if (gouy_vs > 0 && gouy_vs < 1e-8)
    fprintf('  PASS: gouyPhase very small z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gouyPhase very small z\n');
    failed = failed + 1;
end

% testRayleighDistanceVerySmallW0
w0_vs = 1e-9;
zr_vs = PhysicalConstants.rayleighDistance(w0_vs, lambda);
if (zr_vs > 0)
    fprintf('  PASS: rayleighDistance very small w0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: rayleighDistance very small w0\n');
    failed = failed + 1;
end

% testRayleighDistanceVeryLargeW0
w0_vl = 1e-3;
zr_vl = PhysicalConstants.rayleighDistance(w0_vl, lambda);
if (zr_vl > 1)
    fprintf('  PASS: rayleighDistance very large w0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: rayleighDistance very large w0\n');
    failed = failed + 1;
end

% testWaveNumberVerySmallLambda
lambda_vs = 100e-12;
k_vs = PhysicalConstants.waveNumber(lambda_vs);
if (k_vs > 1e10)
    fprintf('  PASS: waveNumber very small lambda\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber very small lambda\n');
    failed = failed + 1;
end

% testWaveNumberVeryLargeLambda
lambda_vl = 1;
k_vl = PhysicalConstants.waveNumber(lambda_vl);
if (k_vl > 0 && k_vl < 10)
    fprintf('  PASS: waveNumber very large lambda\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber very large lambda\n');
    failed = failed + 1;
end

% testPlanckReducedRelation
if (abs(PhysicalConstants.planck / (2*pi) - PhysicalConstants.planck_reduced) < 1e-40)
    fprintf('  PASS: planck reduced relation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: planck reduced relation\n');
    failed = failed + 1;
end

% testImpedanceVacuumRelation
eps0 = PhysicalConstants.vacuum_permittivity;
mu0 = PhysicalConstants.vacuum_permeability;
eta_calc = sqrt(mu0 / eps0);
if (abs(eta_calc - PhysicalConstants.impedance_vacuum) < 1)
    fprintf('  PASS: impedance vacuum relation\n');
    passed = passed + 1;
else
    fprintf('  FAIL: impedance vacuum relation\n');
    failed = failed + 1;
end

% testRadiusOfCurvatureLargeZ
R_l = PhysicalConstants.radiusOfCurvature(10*zr, zr);
if (R_l > zr)
    fprintf('  PASS: radiusOfCurvature large z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: radiusOfCurvature large z\n');
    failed = failed + 1;
end

% testWaistAtZNegativeZ
w_nz = PhysicalConstants.waistAtZ(w0, -0.05, lambda);
if (w_nz > w0)
    fprintf('  PASS: waistAtZ negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waistAtZ negative z\n');
    failed = failed + 1;
end

% testWaveNumberMatrixInput
lambda_mat = [532e-9 633e-9; 1064e-9 1550e-9];
k_mat = PhysicalConstants.waveNumber(lambda_mat);
if (all(all(k_mat > 0)) && isequal(size(k_mat), size(lambda_mat)))
    fprintf('  PASS: waveNumber matrix input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: waveNumber matrix\n');
    failed = failed + 1;
end

fprintf('\n=== PhysicalConstants: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
