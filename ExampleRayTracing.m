%% Canonical Example: Ray Tracing Integration
%% This script demonstrates ray tracing propagation via phase gradients.
% Example script for Ray Tracing Integration
addpath('ParaxialBeams');

% 1. Create a Gaussian Beam
lambda = 632.8e-9;
w0 = 100e-6;
% NEW API: Use Parameters object
params = GaussianParameters(0, w0, lambda);
beam = GaussianBeam(0, params);
zr = params.RayleighDistance;

% 2. Create a Bundle of Rays (Concentric pattern)
Nr = 5; Ntheta = 12; maxR = 2*w0;
bundle = RayBundle.createConcentric(Nr, Ntheta, maxR);

% 3. Propagate the bundle
z_final = 2*zr;
dz = zr/20;
RayTracer.propagate(bundle, beam, z_final, dz, 'RK4');

% 4. Visualize
VisualizationUtils.plotRays3D(bundle, 'b');
VisualizationUtils.plotRays2D(bundle, 'xz', 'r');
