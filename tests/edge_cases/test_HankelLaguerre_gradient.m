% Compatible with GNU Octave and MATLAB
% Regression tests for HankelLaguerre gradient fix.
% Tests the polar-coordinate phase gradient for vortex beams.

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
setpaths();

fprintf('=== HankelLaguerre Gradient Fix Tests ===\n\n');
passed = 0;
failed = 0;

%% Parameters
w0 = 100e-6;
lambda = 632.8e-9;
zr = pi * w0^2 / lambda;

%% ===================================================================
% TEST 1: beamHasVortex detection
%% ===================================================================

% testBeamHasVortexGaussian
beam_g = GaussianBeam(w0, lambda);
if ~RayTracer.beamHasVortex(beam_g)
    fprintf('  PASS: GaussianBeam has no vortex\n');
    passed = passed + 1;
else
    fprintf('  FAIL: GaussianBeam should not have vortex\n');
    failed = failed + 1;
end

% testBeamHasVortexLaguerreL0
beam_l0 = HankelLaguerre(w0, lambda, 0, 0, 1);
if ~RayTracer.beamHasVortex(beam_l0)
    fprintf('  PASS: HankelLaguerre(l=0) has no vortex\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelLaguerre(l=0) should not have vortex\n');
    failed = failed + 1;
end

% testBeamHasVortexLaguerreL8
beam_l8 = HankelLaguerre(w0, lambda, 8, 0, 1);
if RayTracer.beamHasVortex(beam_l8)
    fprintf('  PASS: HankelLaguerre(l=8) has vortex\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelLaguerre(l=8) should have vortex\n');
    failed = failed + 1;
end

% testBeamHasVortexHermite
beam_h = HankelHermite(w0, lambda, 1, 1, 11);
if ~RayTracer.beamHasVortex(beam_h)
    fprintf('  PASS: HankelHermite has no vortex (Cartesian)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelHermite should not have vortex\n');
    failed = failed + 1;
end

%% ===================================================================
% TEST 2: Gradient at vortex-singular angles (θ=0, π)
%% ===================================================================
% The original bug: at θ=0 and θ=π, the Cartesian complex gradient
% produced explosive, asymmetric slopes (e.g., sy ~ 33,000 rad/m).
% With the polar gradient fix, slopes should be FINITE and PHYSICAL.

R_obs = 0.6 * w0;

% testGradientAtTheta0
% At θ=0 (x=R_obs, y=0), the vortex singularity causes Cartesian
% gradient to fail. Polar gradient should be finite.
theta_0 = 0;
x_0 = R_obs * cos(theta_0);
y_0 = R_obs * sin(theta_0);
[sx_0, sy_0] = HankelRayTracer.calculateSlopes(beam_l8, x_0, y_0, 0, 1);
if isfinite(sx_0) && isfinite(sy_0) && abs(sx_0) < 1000 && abs(sy_0) < 1000
    fprintf('  PASS: Gradient at θ=0 is finite (sx=%.2f, sy=%.2f rad/m)\n', sx_0, sy_0);
    passed = passed + 1;
else
    fprintf('  FAIL: Gradient at θ=0 exploded or non-finite (sx=%.2f, sy=%.2f rad/m)\n', sx_0, sy_0);
    failed = failed + 1;
end

% testGradientAtThetaPi
% At θ=π (x=-R_obs, y=0), same singularity issue.
theta_pi = pi;
x_pi = R_obs * cos(theta_pi);
y_pi = R_obs * sin(theta_pi);
[sx_pi, sy_pi] = HankelRayTracer.calculateSlopes(beam_l8, x_pi, y_pi, 0, 1);
if isfinite(sx_pi) && isfinite(sy_pi) && abs(sx_pi) < 1000 && abs(sy_pi) < 1000
    fprintf('  PASS: Gradient at θ=π is finite (sx=%.2f, sy=%.2f rad/m)\n', sx_pi, sy_pi);
    passed = passed + 1;
else
    fprintf('  FAIL: Gradient at θ=π exploded or non-finite (sx=%.2f, sy=%.2f rad/m)\n', sx_pi, sy_pi);
    failed = failed + 1;
end

% testGradientAtThetaPiOver4
% At θ=π/4, the original Cartesian method worked fine.
% Verify the polar gradient also works (no regression).
theta_pi4 = pi/4;
x_pi4 = R_obs * cos(theta_pi4);
y_pi4 = R_obs * sin(theta_pi4);
[sx_pi4, sy_pi4] = HankelRayTracer.calculateSlopes(beam_l8, x_pi4, y_pi4, 0, 1);
if isfinite(sx_pi4) && isfinite(sy_pi4) && abs(sx_pi4) < 1000 && abs(sy_pi4) < 1000
    fprintf('  PASS: Gradient at θ=π/4 is finite (sx=%.2f, sy=%.2f rad/m)\n', sx_pi4, sy_pi4);
    passed = passed + 1;
else
    fprintf('  FAIL: Gradient at θ=π/4 exploded or non-finite (sx=%.2f, sy=%.2f rad/m)\n', sx_pi4, sy_pi4);
    failed = failed + 1;
end

%% ===================================================================
% TEST 3: Symmetry check — gradients at symmetric angles should match
%% ===================================================================
% θ=0 and θ=π should have similar |sy| (mirror symmetry)
% since the beam profile is symmetric.

if abs(abs(sy_0) - abs(sy_pi)) < 100
    fprintf('  PASS: Symmetric gradient at θ=0 and θ=π (|sy|=%.2f vs %.2f)\n', ...
        abs(sy_0), abs(sy_pi));
    passed = passed + 1;
else
    fprintf('  FAIL: Asymmetric gradient at θ=0 and θ=π (|sy|=%.2f vs %.2f)\n', ...
        abs(sy_0), abs(sy_pi));
    failed = failed + 1;
end

%% ===================================================================
% TEST 4: Full propagation — rays should NOT explode
%% ===================================================================
% The original bug caused rays at θ=0 to jump to 25+ meters in one step.
% After fix, rays should propagate normally.

% testPropagationNoExplosion
bundle = RayBundle.createCircularContour(8, R_obs);
bundle.ht(:) = 1;
dz = zr / 64;
z_final = zr / 10;  % 10 steps

bundle = HankelRayTracer.propagate(bundle, beam_l8, z_final, dz, 'RK4');

x_final = bundle.x(:);
y_final = bundle.y(:);
max_displacement = max(sqrt(x_final.^2 + y_final.^2));

% After 10 steps, rays should still be within ~10x initial radius, not meters
if max_displacement < 10 * R_obs && all(isfinite(x_final)) && all(isfinite(y_final))
    fprintf('  PASS: Propagation stable (max r = %.2e m, not 25+ m)\n', max_displacement);
    passed = passed + 1;
else
    fprintf('  FAIL: Propagation exploded (max r = %.2e m)\n', max_displacement);
    failed = failed + 1;
end

% testPropagationDivergence
% Verify rays still diverge as expected (ray bundle should spread)
rms_initial = sqrt(mean(bundle.x(:,:,1).^2 + bundle.y(:,:,1).^2));
rms_final = sqrt(mean(bundle.x(:,:,end).^2 + bundle.y(:,:,end).^2));
if rms_final > rms_initial
    fprintf('  PASS: Rays diverge as expected (RMS: %.2e → %.2e m)\n', rms_initial, rms_final);
    passed = passed + 1;
else
    fprintf('  FAIL: Rays did not diverge (RMS: %.2e → %.2e m)\n', rms_initial, rms_final);
    failed = failed + 1;
end

%% ===================================================================
% TEST 5: l=0 HankelLaguerre (no vortex) — no regression
%% ===================================================================
% l=0 has no vortex, should still work with Cartesian gradient.

bundle_l0 = RayBundle.createCircularContour(4, R_obs);
bundle_l0.ht(:) = 1;
beam_l0_test = HankelLaguerre(w0, lambda, 0, 0, 1);

bundle_l0 = HankelRayTracer.propagate(bundle_l0, beam_l0_test, z_final, dz, 'RK4');

if all(isfinite(bundle_l0.x(:))) && all(isfinite(bundle_l0.y(:)))
    fprintf('  PASS: l=0 HankelLaguerre propagation stable (no regression)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: l=0 HankelLaguerre propagation failed\n');
    failed = failed + 1;
end

%% ===================================================================
% TEST 6: HankelHermite (no vortex) — no regression
%% ===================================================================
beam_h_test = HankelHermite(w0, lambda, 1, 1, 11);
bundle_h = RayBundle.createGrid(3, 3, w0, w0);
bundle_h.ht(:) = 11;

bundle_h = HankelRayTracer.propagate(bundle_h, beam_h_test, z_final, dz, 'RK4');

if all(isfinite(bundle_h.x(:))) && all(isfinite(bundle_h.y(:)))
    fprintf('  PASS: HankelHermite propagation stable (no regression)\n');
    passed = passed + 1;
else
    fprintf('  FAIL: HankelHermite propagation failed\n');
    failed = failed + 1;
end

%% ===================================================================
% SUMMARY
%% ===================================================================
fprintf('\n=== HankelLaguerre Gradient Fix: %d/%d passed ===\n', passed, passed + failed);

if failed ~= 0
    error('Tests failed: %d/%d', failed, passed + failed);
end