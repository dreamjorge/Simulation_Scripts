%% Gaussian Beam Script (Refactored)
% Using new utility classes
% Ugalde-Ontiveros J.A.

addpath ParaxialBeams
addpath ParaxialBeams\Addons
addpath ParaxialBeams\Addons\export_fig-master
addpath ParaxialBeams\Addons\Plots_Functions

%% Physical parameters [meters]
InitialWaist = 100e-6;
Wavelength = 632.8e-9;

%% Using PhysicalConstants
PC = PhysicalConstants;
k = PC.waveNumber(Wavelength);
zr = PC.rayleighDistance(InitialWaist, Wavelength);

fprintf('Physical Parameters:\n');
fprintf('  Wavelength: %g nm\n', Wavelength*1e9);
fprintf('  InitialWaist: %g um\n', InitialWaist*1e6);
fprintf('  Rayleigh distance: %g m\n', zr);

%% Create grid using GridUtils
Nx = 2^10;
Nz = 2^8;
Dz = 2*zr;

grid = GridUtils(Nx, Nx, 1, 1, Nz, Dz);
maxWaist = PC.waistAtZ(InitialWaist, Dz, Wavelength, zr);
grid.Dx = 1.2 * 2 * maxWaist;
grid.dx = grid.Dx / Nx;

[X, Y] = grid.create2DGrid();
[Kx, Ky] = grid.createFreqGrid();

fprintf('\nGrid Parameters:\n');
fprintf('  Nx = %d, Nz = %d\n', Nx, Nz);
fprintf('  Dx = %g m\n', grid.Dx);

%% FFT utilities
fftOps = FFTUtils(true, true);

%% Gaussian Beam
GaussianBeamParameters = GaussianParameters(0, InitialWaist, Wavelength);
GB = GaussianBeam( hypot(X,Y), GaussianBeamParameters);

% Plot
figure(1)
plotOpticalField(X, Y, abs(GB.OpticalField).^2, mapgreen, 'microns');

%% Propagate using angular spectrum
g0 = GB.OpticalField;
zPlanes = linspace(0, Dz, Nz);

for iz = 1:Nz
    z = zPlanes(iz);
    H = fftOps.transferFunction(Kx, Ky, z, Wavelength);
    G = fftOps.fft2(g0);
    gProp = fftOps.ifft2(G .* H);
    % Store or plot gProp
end

fprintf('\nPropagation complete\n');
