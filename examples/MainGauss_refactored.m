%% Canonical Example: Gaussian Beam Propagation
%% This script demonstrates the recommended workflow for beam simulation.
% Gaussian Beam propagation using new utility classes
% Ugalde-Ontiveros J.A.
%
% This version uses:
%   - PhysicalConstants for physical parameter calculations
%   - GridUtils for grid generation
%   - FFTUtils for normalized FFT operations

%% Add paths
addpath(fullfile('ParaxialBeams'))
addpath(fullfile('ParaxialBeams', 'Addons'))
addpath(fullfile('ParaxialBeams', 'Addons', 'export_fig-master'))
addpath(fullfile('ParaxialBeams', 'Addons', 'panel-2.14'))
addpath(fullfile('ParaxialBeams', 'Addons', 'Plots_Functions'))

%% Physical parameters [microns -> meters]
InitialWaist = 100e-6;      % 100 microns
Wavelength  = 632.8e-9;    % HeNe laser (632.8 nm)

%% Calculate derived parameters using PhysicalConstants
PC = PhysicalConstants;
k = PC.waveNumber(Wavelength);
zr = PC.rayleighDistance(InitialWaist, Wavelength);

fprintf('Physical Parameters:\n');
fprintf('  Wavelength: %g nm\n', Wavelength*1e9);
fprintf('  InitialWaist: %g um\n', InitialWaist*1e6);
fprintf('  Rayleigh distance: %g m\n', zr);
fprintf('  Wave number k: %g 1/m\n', k);

%% Create grid using GridUtils
Nx = 2^10;
Dz = 2*zr;
Nz = 2^8;

% Create grid instance
grid = GridUtils(Nx, Nx, 1, 1, Nz, Dz);

% Estimate optimal window size
maxWaist = PC.waistAtZ(InitialWaist, Dz, Wavelength, zr);
grid.Dx = 1.2 * 2 * maxWaist;
grid.Dy = grid.Dx;
grid.dx = grid.Dx / Nx;
grid.dy = grid.Dy / Nx;

% Get grid arrays
[X, Y] = grid.create2DGrid();
[Kx, Ky] = grid.createFreqGrid();

fprintf('\nGrid Parameters:\n');
fprintf('  Nx = %d, Nz = %d\n', Nx, Nz);
fprintf('  Dx = %g m, Dz = %g m\n', grid.Dx, grid.Dz);
fprintf('  dx = %g m\n', grid.dx);

%% FFT utilities
fftOps = FFTUtils(true, true);  % normalize=true, shift=true

% Colormap for visualization
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);

%% Gaussian Beam at z = 0
GaussianBeamParameters = GaussianParameters(0, InitialWaist, Wavelength);
[R, ~] = cart2pol(X, Y);
GB = GaussianBeam(R, GaussianBeamParameters);

% Plot field
figure(1)
plotOpticalField(X, Y, abs(GB.OpticalField).^2, mapgreen, 'microns');
plotCircle(0, 0, GaussianBeamParameters.InitialWaist);

%% Compare old vs new FFT approach
% Old way (scattered in codebase):
%   G = fftshift(fft2(ifftshift(g)));
%   g = fftshift(ifft2(ifftshift(G)));

% New way (consistent):
%   G = fftOps.fft2(g);
%   g = fftOps.ifft2(G);

%% Example: Propagate beam using angular spectrum
% Create initial field
g0 = GB.OpticalField;

% Propagate to multiple z planes
zPlanes = linspace(0, Dz, Nz);
fields = zeros(Nx, Nx, Nz);

for iz = 1:Nz
    z = zPlanes(iz);
    
    % Compute transfer function
    H = fftOps.transferFunction(Kx, Ky, z, Wavelength);
    
    % Propagate
    G = fftOps.fft2(g0);
    gProp = fftOps.ifft2(G .* H);
    fields(:, :, iz) = gProp;
end

fprintf('\nPropagation complete: %d z-planes computed\n', Nz);

%% Export parameters info
info = struct('InitialWaist', InitialWaist, ...
              'Wavelength', Wavelength, ...
              'RayleighDistance', zr, ...
              'Nx', Nx, 'Nz', Nz, ...
              'Dx', grid.Dx, 'Dz', Dz);
          
save('GaussianSimulationInfo.mat', 'info');
fprintf('\nInfo saved to GaussianSimulationInfo.mat\n');
