%% Canonical Example: Multi-Mode Beam Demonstration
%% This script demonstrates multi-mode Hermite and Laguerre beams.
% Using refactored class architecture
% Ugalde-Ontiveros J.A.

addpath(fullfile('ParaxialBeams'))

%% Parameters
w0 = 100e-6;
lambda = 632.8e-9;
z = 0.5; % Propagate to 0.5m

%% Grid
Nx = 512;
Dx = 1.5e-3;
grid = GridUtils(Nx, Nx, Dx, Dx);
[X, Y] = grid.create2DGrid();
[R, Theta] = cart2pol(X, Y);

%% Hermite-Gaussian (1,1)
fprintf('Calculating Hermite-Gaussian (1,1)...\n');
hp = HermiteParameters(z, w0, lambda, 1, 1);
hb = HermiteBeam(X, Y, hp);

%% Laguerre-Gaussian (0,1)
fprintf('Calculating Laguerre-Gaussian (0,1)...\n');
lp = LaguerreParameters(z, w0, lambda, 1, 0); % l=1, p=0
lb = LaguerreBeam(R, Theta, lp);

%% Visualization
figure(1);
subplot(1,2,1);
imagesc(abs(hb.OpticalField).^2);
colormap('hot'); axis square; title('HG (1,1) Intensity');

subplot(1,2,2);
imagesc(abs(lb.OpticalField).^2);
colormap('hot'); axis square; title('LG (0,1) Intensity');

fprintf('Done.\n');
