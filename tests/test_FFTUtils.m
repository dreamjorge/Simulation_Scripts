#!/usr/bin/env octave
% Tests for FFTUtils

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'ParaxialBeams'));

fprintf('=== FFTUtils Tests ===\n\n');
passed = 0;
failed = 0;

lambda = 632.8e-9;

% Setup
[X, Y] = meshgrid(linspace(-1,1,64), linspace(-1,1,64));
R = sqrt(X.^2 + Y.^2);
g = exp(-R.^2);
fftOps = FFTUtils(true, true);

% test fft2 roundtrip
g_rec = fftOps.ifft2(fftOps.fft2(g));
if (max(max(abs(g - g_rec))) < 1e-10)
    fprintf('  PASS: fft2 roundtrip\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2 roundtrip\n');
    failed = failed + 1;
end

% test transferFunction at z=0
[Kx, Ky] = meshgrid(linspace(-1e6,1e6,32));
H = FFTUtils.transferFunction(Kx, Ky, 0, lambda);
if (max(max(abs(H - 1))) < 1e-10)
    fprintf('  PASS: transferFunction at z=0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferFunction at z=0\n');
    failed = failed + 1;
end

% test fftn
g3d = exp(-(X.^2 + Y.^2));
g3d_rec = fftOps.ifftn(fftOps.fftn(g3d));
if (max(max(max(abs(g3d - g3d_rec)))) < 1e-10)
    fprintf('  PASS: fftn 3D roundtrip\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fftn 3D roundtrip\n');
    failed = failed + 1;
end

% test propagate
g_prop = fftOps.propagate(g, zeros(64), zeros(64), 0, lambda);
if (max(max(abs(g - g_prop))) < 1e-10)
    fprintf('  PASS: propagate roundtrip\n');
    passed = passed + 1;
else
    fprintf('  FAIL: propagate roundtrip\n');
    failed = failed + 1;
end

% test fft2_centered
G_c = FFTUtils.fft2_centered(g);
g_c_rec = FFTUtils.ifft2_centered(G_c);
if (max(max(abs(g - g_c_rec))) < 1e-10)
    fprintf('  PASS: fft2_centered\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2_centered\n');
    failed = failed + 1;
end

% test transferFunction unit magnitude
H_mag = FFTUtils.transferFunction(0, 0, 0.1, lambda);
if (abs(abs(H_mag) - 1) < 1e-10)
    fprintf('  PASS: transferFunction unit magnitude\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferFunction magnitude\n');
    failed = failed + 1;
end

fprintf('\n=== FFTUtils: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
