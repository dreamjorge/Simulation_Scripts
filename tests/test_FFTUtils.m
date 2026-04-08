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

% test transferSimple
H_s = FFTUtils.transferSimple(Kx, Ky, 0.1, lambda);
if (numel(H_s) == numel(Kx))
    fprintf('  PASS: transferSimple\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferSimple\n');
    failed = failed + 1;
end

% test fft2 without shift
fft_noshift = FFTUtils(true, false);
g_noshift = fft_noshift.fft2(g);
if (numel(g_noshift) == numel(g))
    fprintf('  PASS: fft2 without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2 no shift\n');
    failed = failed + 1;
end

% test ifft2 without shift
g_rec_noshift = fft_noshift.ifft2(g_noshift);
if (max(max(abs(g - g_rec_noshift))) < 1e-6)
    fprintf('  PASS: ifft2 without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ifft2 no shift\n');
    failed = failed + 1;
end

% test fftn without shift
fft3d_noshift = FFTUtils(true, false);
g3d_noshift = fft3d_noshift.fftn(g3d);
if (numel(g3d_noshift) == numel(g3d))
    fprintf('  PASS: fftn without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fftn no shift\n');
    failed = failed + 1;
end

% test ifftn without shift
g3d_rec_noshift = fft3d_noshift.ifftn(g3d_noshift);
if (max(max(max(abs(g3d - g3d_rec_noshift)))) < 1e-6)
    fprintf('  PASS: ifftn without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ifftn no shift\n');
    failed = failed + 1;
end

% test normalize false
fft_nonorm = FFTUtils(false, true);
g_nonorm = fft_nonorm.fft2(g);
g_rec_nonorm = fft_nonorm.ifft2(g_nonorm);
if (max(max(abs(g - g_rec_nonorm))) < 1e-10)
    fprintf('  PASS: normalize false\n');
    passed = passed + 1;
else
    fprintf('  FAIL: normalize false\n');
    failed = failed + 1;
end

% test fft2_std
G_std = FFTUtils.fft2_std(g);
g_std_rec = FFTUtils.ifft2_std(G_std);
if (max(max(abs(g - g_std_rec))) < 1e-10)
    fprintf('  PASS: fft2_std\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2_std\n');
    failed = failed + 1;
end

% test propagate with kx ky
[Kx_test, Ky_test] = meshgrid(linspace(-1e5,1e5,16), linspace(-1e5,1e5,16));
g_prop_k = fftOps.propagate(g(1:16,1:16), Kx_test, Ky_test, 0.01, lambda);
if (numel(g_prop_k) == 256)
    fprintf('  PASS: propagate with kx ky\n');
    passed = passed + 1;
else
    fprintf('  FAIL: propagate kx ky\n');
    failed = failed + 1;
end

% test transferFunction complex
H_complex = FFTUtils.transferFunction(Kx(1:8,1:8), Ky(1:8,1:8), 0.05, lambda);
if (~isreal(H_complex))
    fprintf('  PASS: transferFunction complex\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferFunction complex\n');
    failed = failed + 1;
end

% test transferSimple at origin
H_s_0 = FFTUtils.transferSimple(0, 0, 0.1, lambda);
if (abs(H_s_0 - exp(1i * 2*pi/lambda * 0.1)) < 1e-10)
    fprintf('  PASS: transferSimple at origin\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferSimple origin\n');
    failed = failed + 1;
end

% test fft2_centered output size
G_cc = FFTUtils.fft2_centered(g);
if (size(G_cc) == size(g))
    fprintf('  PASS: fft2_centered size\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2_centered size\n');
    failed = failed + 1;
end

% test ifft2_centered roundtrip
g_cc_rt = FFTUtils.ifft2_centered(G_cc);
if (max(max(abs(g - g_cc_rt))) < 1e-10)
    fprintf('  PASS: ifft2_centered roundtrip\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ifft2_centered rt\n');
    failed = failed + 1;
end

% test propagate z>0
g_prop_z = fftOps.propagate(g, zeros(64), zeros(64), 0.1, lambda);
if (all(all(isfinite(g_prop_z))))
    fprintf('  PASS: propagate z>0\n');
    passed = passed + 1;
else
    fprintf('  FAIL: propagate z>0\n');
    failed = failed + 1;
end

% test FFTUtils constructor default
fft_def = FFTUtils();
if (fft_def.normalize == true && fft_def.shiftFlag == true)
    fprintf('  PASS: constructor default\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor default\n');
    failed = failed + 1;
end

% test FFTUtils constructor normalize only
fft_norm = FFTUtils(false);
if (fft_norm.normalize == false && fft_norm.shiftFlag == true)
    fprintf('  PASS: constructor normalize only\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor normalize\n');
    failed = failed + 1;
end

% test propagate negative z
g_prop_nz = fftOps.propagate(g, zeros(64), zeros(64), -0.01, lambda);
if (all(all(isfinite(g_prop_nz))))
    fprintf('  PASS: propagate negative z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: propagate negative z\n');
    failed = failed + 1;
end

% test transferFunctionLargeZ
[H_lz] = FFTUtils.transferFunction(Kx, Ky, 1, lambda);
if (numel(H_lz) == numel(Kx))
    fprintf('  PASS: transferFunction large z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferFunction large z\n');
    failed = failed + 1;
end

% test transferSimpleLargeKx
Kx_l = linspace(-1e7, 1e7, 32);
Ky_l = zeros(1, 32);
H_sl = FFTUtils.transferSimple(Kx_l, Ky_l, 0.1, lambda);
if (numel(H_sl) == 32)
    fprintf('  PASS: transferSimple large kx\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferSimple large kx\n');
    failed = failed + 1;
end

% test fftn roundtrip
g_3d = rand(8, 8, 8);
fft3d = FFTUtils(true, true);
g_3d_rec = fft3d.ifftn(fft3d.fftn(g_3d));
if (max(max(max(abs(g_3d - g_3d_rec)))) < 1e-10)
    fprintf('  PASS: fftn roundtrip\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fftn roundtrip\n');
    failed = failed + 1;
end

% test fftn without shift
fft3d_ns = FFTUtils(true, false);
g_3d_ns = fft3d_ns.fftn(g_3d);
if (numel(g_3d_ns) == numel(g_3d))
    fprintf('  PASS: fftn without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fftn no shift\n');
    failed = failed + 1;
end

% test ifftn without shift
g_3d_ns_rec = fft3d_ns.ifftn(g_3d_ns);
if (max(max(max(abs(g_3d - g_3d_ns_rec)))) < 1e-6)
    fprintf('  PASS: ifftn without shift\n');
    passed = passed + 1;
else
    fprintf('  FAIL: ifftn no shift\n');
    failed = failed + 1;
end

% test fft2 zero input
g_z = zeros(16, 16);
g_z_rec = fftOps.ifft2(fftOps.fft2(g_z));
if (max(max(abs(g_z - g_z_rec))) < 1e-10)
    fprintf('  PASS: fft2 zero input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2 zero input\n');
    failed = failed + 1;
end

% test fft2 ones input
g_o = ones(16, 16);
g_o_rec = fftOps.ifft2(fftOps.fft2(g_o));
if (max(max(abs(g_o - g_o_rec))) < 1e-10)
    fprintf('  PASS: fft2 ones input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2 ones input\n');
    failed = failed + 1;
end

% test propagate with small z
g_prop_sz = fftOps.propagate(g, zeros(64), zeros(64), 1e-6, lambda);
if (all(all(isfinite(g_prop_sz))))
    fprintf('  PASS: propagate small z\n');
    passed = passed + 1;
else
    fprintf('  FAIL: propagate small z\n');
    failed = failed + 1;
end

% test transferFunctionZeroKxKy
H_zk = FFTUtils.transferFunction(zeros(8), zeros(8), 0.1, lambda);
if (abs(H_zk - exp(1i * 2*pi/lambda * 0.1)) < 1e-10)
    fprintf('  PASS: transferFunction zero kx ky\n');
    passed = passed + 1;
else
    fprintf('  FAIL: transferFunction zero kx ky\n');
    failed = failed + 1;
end

% test fft2_centered with complex input
g_cx = exp(1i*rand(16,16));
G_cx = FFTUtils.fft2_centered(g_cx);
g_cx_rec = FFTUtils.ifft2_centered(G_cx);
if (max(max(abs(g_cx - g_cx_rec))) < 1e-10)
    fprintf('  PASS: fft2_centered complex input\n');
    passed = passed + 1;
else
    fprintf('  FAIL: fft2_centered complex\n');
    failed = failed + 1;
end

% test FFTUtils constructor both false
fft_ff = FFTUtils(false, false);
if (fft_ff.normalize == false && fft_ff.shiftFlag == false)
    fprintf('  PASS: constructor both false\n');
    passed = passed + 1;
else
    fprintf('  FAIL: constructor both false\n');
    failed = failed + 1;
end

fprintf('\n=== FFTUtils: %d/%d passed ===\n', passed, passed + failed);

if (failed == 0)
    exit(0);
else
    exit(1);
end
