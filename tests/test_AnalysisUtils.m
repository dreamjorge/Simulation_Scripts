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

fprintf('\n=== AnalysisUtils: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
