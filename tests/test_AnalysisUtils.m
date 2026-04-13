#!/usr/bin/env octave
% Tests for AnalysisUtils

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== AnalysisUtils Tests ===\n\n');
passed = 0;
failed = 0;

% testGradientRZ
fr = ones(1, 10); fz = ones(1, 10);
mzr = AnalysisUtils.gradientRZ(fr, fz, 1e7, 1e-4, 1e-4, 0, 0);
if (isfinite(mzr))
    fprintf('  PASS: gradientRZ\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ\n');
    failed = failed + 1;
end

% testGradientRZNonUniform
fr_n = exp(-((1:10)*1e-4).^2);
fz_n = exp(-((1:10)*1e-4).^2);
mzr_n = AnalysisUtils.gradientRZ(fr_n, fz_n, 1e7, 1e-4, 1e-4, 5e-5, 5e-5);
if (isfinite(mzr_n))
    fprintf('  PASS: gradientRZ non-uniform\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ non-uniform\n');
    failed = failed + 1;
end

% testGradientXYZ
fyz = ones(16, 16); fxz = ones(16, 16); fxy = ones(16, 16);
[mzx, mzy, mxy] = AnalysisUtils.gradientXYZ(fyz, fxz, fxy, 1e6, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5);
if (numel(mzx) > 0)
    fprintf('  PASS: gradientXYZ\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ\n');
    failed = failed + 1;
end

% testGradientXYZGaussian
[Xg, Yg] = meshgrid(linspace(-1,1,16), linspace(-1,1,16));
fyz_g = exp(-(Xg.^2 + Yg.^2));
fxz_g = exp(-(Xg.^2 + Yg.^2));
fxy_g = exp(-(Xg.^2 + Yg.^2));
[mzx_g, mzy_g, mxy_g] = AnalysisUtils.gradientXYZ(fyz_g, fxz_g, fxy_g, 1e7, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5);
if (numel(mzx_g) > 0)
    fprintf('  PASS: gradientXYZ with Gaussian\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ Gaussian\n');
    failed = failed + 1;
end

% testGradientXYZBounds
fyz_b = ones(32, 32); fxz_b = ones(32, 32); fxy_b = ones(32, 32);
try
    [mzx_b, mzy_b, mxy_b] = AnalysisUtils.gradientXYZ(fyz_b, fxz_b, fxy_b, 1e6, 1e-4, 1e-4, 1e-4, 5e-5, 5e-5, 5e-5);
    fprintf('  PASS: gradientXYZ bounds\n');
    passed = passed + 1;
catch
    fprintf('  FAIL: gradientXYZ bounds\n');
    failed = failed + 1;
end

% testGradientRZAtSpecificPoint
fr_s = linspace(0, 1, 20);
fz_s = linspace(0, 1, 20);
mzr_s = AnalysisUtils.gradientRZ(fr_s, fz_s, 1e7, 1e-4, 1e-4, 0.0005, 0.0005);
if (isfinite(mzr_s))
    fprintf('  PASS: gradientRZ at specific point\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ at point\n');
    failed = failed + 1;
end

% testGradientXYZAtSpecificPoint
fyz_p = ones(20, 20); fxz_p = ones(20, 20); fxy_p = ones(20, 20);
[mzx_p, mzy_p, mxy_p] = AnalysisUtils.gradientXYZ(fyz_p, fxz_p, fxy_p, 1e6, 1e-4, 1e-4, 1e-4, 0.0005, 0.0005, 0.0005);
if (isfinite(mzx_p) && isfinite(mzy_p))
    fprintf('  PASS: gradientXYZ at specific point\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ at point\n');
    failed = failed + 1;
end

% testGradientRZWithZeroField
fr_z = zeros(1, 10);
fz_z = ones(1, 10);
try
    mzr_z = AnalysisUtils.gradientRZ(fr_z, fz_z, 1e7, 1e-4, 1e-4, 0, 0);
    fprintf('  PASS: gradientRZ with zero field\n');
    passed = passed + 1;
catch
    fprintf('  FAIL: gradientRZ zero field\n');
    failed = failed + 1;
end

% testGradientXYZNonUniform
fyz_n = rand(16, 16);
fxz_n = rand(16, 16);
fxy_n = rand(16, 16);
[mzx_n, mzy_n, mxy_n] = AnalysisUtils.gradientXYZ(fyz_n, fxz_n, fxy_n, 1e6, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5);
if (numel(mzx_n) > 0)
    fprintf('  PASS: gradientXYZ non-uniform\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ non-uniform\n');
    failed = failed + 1;
end

% testGradientXYZOutputSize
[mzx_s, mzy_s, mxy_s] = AnalysisUtils.gradientXYZ(ones(8,8), ones(8,8), ones(8,8), 1e6, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5);
if (numel(mzx_s) == 1)
    fprintf('  PASS: gradientXYZ scalar output\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ scalar\n');
    failed = failed + 1;
end

% testGradientRZLargeK
mzr_lk = AnalysisUtils.gradientRZ(ones(1,10), ones(1,10), 1e9, 1e-4, 1e-4, 0, 0);
if (isfinite(mzr_lk))
    fprintf('  PASS: gradientRZ large k\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ large k\n');
    failed = failed + 1;
end

% testGradientXYZLargeK
[mzx_lk, mzy_lk, mxy_lk] = AnalysisUtils.gradientXYZ(ones(16,16), ones(16,16), ones(16,16), 1e9, 1e-4, 1e-4, 1e-4, 1e-5, 1e-5, 1e-5);
if (numel(mzx_lk) > 0)
    fprintf('  PASS: gradientXYZ large k\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ large k\n');
    failed = failed + 1;
end

% testGradientRZSmallStep
mzr_ss = AnalysisUtils.gradientRZ(ones(1,10), ones(1,10), 1e7, 1e-6, 1e-6, 1e-7, 1e-7);
if (isfinite(mzr_ss))
    fprintf('  PASS: gradientRZ small step\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ small step\n');
    failed = failed + 1;
end

% testGradientXYZSmallStep
[mzx_ss, mzy_ss, mxy_ss] = AnalysisUtils.gradientXYZ(ones(16,16), ones(16,16), ones(16,16), 1e6, 1e-6, 1e-6, 1e-6, 1e-7, 1e-7, 1e-7);
if (numel(mzx_ss) > 0)
    fprintf('  PASS: gradientXYZ small step\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientXYZ small step\n');
    failed = failed + 1;
end

% testGradientRZEdgeCase
fr_e = linspace(0, 10, 20);
fz_e = linspace(0, 10, 20);
mzr_e = AnalysisUtils.gradientRZ(fr_e, fz_e, 1e7, 1e-4, 1e-4, 0.0009, 0.0009);
if (isfinite(mzr_e))
    fprintf('  PASS: gradientRZ edge case\n');
    passed = passed + 1;
else
    fprintf('  FAIL: gradientRZ edge\n');
    failed = failed + 1;
end

fprintf('\n=== AnalysisUtils: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end
