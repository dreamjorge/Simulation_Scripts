% app/main_ray_tracing.m
clear;
clc;

% Create a Hermite-Gaussian beam
beam = HermiteGaussianBeam(633e-9, 1e-3, 0, 1, 1);

% Define propagation distances
z_values = linspace(0, 2, 10);  % From z = 0 to z = 2 meters
grid_size = [-5, 5];  % Transverse plane grid

% Perform ray tracing
ray_tracing(beam, z_values, grid_size);
