
%% add path for classes and functions
addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam

%-------------------- indices of Laguerre Gaussian Beams -----------------%
nu      = 6;
mu      = 7;

% Physical parameters [microns]
HermiteInitialWaist  = 10;
InitialWaist         = HermiteInitialWaist/sqrt(nu+mu+1);
Wavelength           = 0.6328;

% Calculating gaussian parameters
GP                   = GaussianParameters(0,InitialWaist,Wavelength);
k                    = GP.k;
RayleighDistance     = GP.RayleighDistance;

% normalized parameters
% 
% Wavelength           = pi;
% InitialWaist         = 1;
% HermiteInitialWaist  = InitialWaist*sqrt(nu+mu+1);
% GP                   = GaussianParameters(0,InitialWaist,Wavelength);
% k                    = GP.k;
% RayleighDistance     = GP.RayleighDistance;
% %------------------------ sampling of vectors ----------------------------%
%First we estimate samplig in z-direction with propagation distance 
% z-direction
Dz     = 2*RayleighDistance;        % z-window (propagation distance)
Nz     = 2^8;                       % number of points in z-direction
dz     = Dz/Nz;                     % Resolution in z
z      = 0:dz:Dz;                   % z-vector z of propagation 

% waist of Laguerre Gauss Beam until z-propagation
maxIndex = max(nu,mu);
MaxHermiteWaist = HermiteBeam.waistHermite(z(end),InitialWaist,RayleighDistance,maxIndex);

%Second we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
N      =  2^9;                   % Number of points in x,y axis
n      = -N/2+.05:N/2-1+.05;      % vector with N-points with resolution 1
Dx     = 2*MaxHermiteWaist;      % Size of window 
dx     = Dx/N;                    % Resolution
x      = n*dx;                    % Vector
y      = x;
[X,Y]  = meshgrid(x,y);

%Last we estimate vectors of frequency for Fourier Transforms
Du     = 1/dx;                    % Size of window 
du     = 1/Dx;                    % Resolution
u      = n*du;                    % Vector
[U]    = meshgrid(u);
%kx,ky vectors
kx     = 2*pi*u;
[Kx]   = meshgrid(kx);

%% ----------------------- Laguerre Gauss in z = 0 --------------------- %%
HG  = HermiteBeam(X,Y,0,InitialWaist,Wavelength,nu,mu);
% Optic Field to propagate 
g   = HG.OpticalField;
% Plot of Function
figure(1)
pcolor(x/(sqrt(2)*InitialWaist), x/(sqrt(2)*InitialWaist), abs(g).^2)
axis square
shading flat
colormap(mapgreen)
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
% Max Peak
pxy = max(max(g));

%% ----------------- Obstruction on Lagurre in z = 0 ------------------- %%
lo      = HermiteInitialWaist/8;    % size of obstruction in terms of waist of Laguerre
xt      = 0;                        % traslation of obstruction in x-axis
yt      = 0;                        % traslation of onstruction in y-axis
[~,rho] = cart2pol(X-xt,X'-yt);     % Convert this in polar coordinates
obo     = double(rho<=lo);          % Create Obstruction   
clear rho                           % Clean Matrix of polar coordinates
% Applying obstruction in optic field
g       = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(2)
pcolor(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)



%%  ----------------------- Physical Propagation ------------------------ %
% paraxial propagator 
prop = exp(1i*dz*(Kx.^2+(Kx').^2)/(2*k));
figure(5)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(N,length(z)); 
gy      = zeros(N,length(z));
% Save field in z = 0 
gx(:,1) = g(N/2+1,:);
gy(:,1) = g(:,N/2+1);


for ii = 2:length(z) % propagation with respect to z
    
    % field before propagation i.e in z(ii-1)
    pxyz         = g(1,1);
    g(1,1)       = pxy;
    fig          = figure(6);

    fig.Position = [ 239 135 1354 733];

%     set(gca,'un','n','pos',[0,0,1,1])
    imagesc(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    title(['z = ', num2str(z(ii))])
    drawnow 
  
   pause(.25)


    %------------------------ End calculating rays -----------------------%   
    %propagating field
    % it's needed correction in phase of FFT
    G = fftshift(fft2((g)));...*exp(-1j.*(u(1)).*(Kx)).*exp(-1j.*(u(1)).*(Kx'));
    %obtain new propagated field
    g = (ifft2(ifftshift(G.*prop)));...*exp(-1j.*(u(1)).*(X))*exp(-1j.*(u(1)).*(X'));
    figure(7)
    imagesc(angle(G))
    %G = G.*exp(1i*pi*50);
    figure(8)
    imagesc(angle(G))
    %obtain new propagated field
    %saving transversal fields
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);
%     
end

