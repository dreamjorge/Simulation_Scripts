% app/main.m
clear;
clc;

% Define beam parameters
wavelength = 633e-9;  % Wavelength in meters
waist = 1e-3;         % Waist size in meters
position = 0;         % Starting position (z = 0)

% Create a Gaussian beam
beam = ParaxialBeam(wavelength, waist, position);

% Pass the beam to the BeamService
beamService = BeamService(beam);

% Calculate properties at z = 1 meter
z = 1;
waist_at_z = beamService.calculateBeamWaist(z);
fprintf('Beam waist at z = %.2f m: %.2e m\n', z, waist_at_z);

% Visualize intensity at z = 1 meter
grid_size = [-5, 5];
plotBeamIntensity(beam, z, grid_size);
