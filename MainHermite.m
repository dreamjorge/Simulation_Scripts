%%                      Main Hermite
%% add path for classes and functions
clear all
addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam

%% indices of Hermite Gaussian Beams 
nu      = 12;
mu      = 11;

%% Physical parameters [microns]
InitialWaist          = 100;
Wavelength            = 0.6328;

HPz0                  = HermiteParameters(0,InitialWaist,Wavelength,nu,mu);

k                     = HPz0.k;
RayleighDistance      = HPz0.RayleighDistance;

HermiteInitialWaistX  = HPz0.HermiteWaistX;
HermiteInitialWaistY  = HPz0.HermiteWaistY;
HermiteInitialWaist   = HPz0.HermiteWaist;

%% normalized parameters
% 
% Wavelength           = pi;
% InitialWaist         = 1;
% HermiteInitialWaist  = InitialWaist*sqrt(nu+mu+1);
% GP                   = GaussianParameters(0,InitialWaist,Wavelength);
% k                    = GP.k;
% RayleighDistance     = GP.RayleighDistance;

%% sampling of vectors 
%First we estimate samplig in z-direction with propagation distance 
% z-direction
Dz    = RayleighDistance/2;     % z-window (propagation distance)
Nz    = 2^6;                  % number of points in z-direction
dz    = Dz/Nz;                % Resolution in z
nz    = 0:Nz-1;               % vector with N-points with resolution 1
z     = nz*dz;                % z-vector z of propagation 

% waist of Hermite Gauss Beam until z-propagation
MaxHermiteWaist = HermiteParameters.getWaist(z(end),InitialWaist,RayleighDistance,nu,mu);

%Second we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
Nx    =  2^8;                % Number of points in x,y axis
n     = -Nx/2:Nx/2-1;         % vector with N-points with resolution 1
Dx    = MaxHermiteWaist;      % Size of window 
dx    = Dx/Nx;                % Resolution
x     = n*dx;                 % Vector
y     = x;
[X,Y] = meshgrid(x,y);

%Last we estimate vectors of frequency for Fourier Transforms
Du    = 1/dx;                 % Size of window 
du    = 1/Dx;                 % Resolution
u     = n*du;                 % Vector
[U]   = meshgrid(u);
%kx,ky vectors
kx    = 2*pi*u;
[Kx]  = meshgrid(kx);

%diferential
dr    = [dx,dx,dz];
 
%% Hermite Gauss in z = 0
HGB   = HermiteBeam(X,Y,HPz0);

% copy of parameters
HPz   = copy(HPz0);

HPz.zCoordinate = z;  

% Optic Field to propagate 
g     = HGB.OpticalFieldHermite;
% Plot of Function
figure(1)
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');

%% Obstruction on Hermite in z = 0

% Obstruction
lx    = HPz0.HermiteWaist/5;
ly    = lx;
xt    = 0;...HPz0.HermiteWaist/5;
yt    = xt;
obx   = double(abs(x-xt)<=lx/2);
oby   = double(abs(x-yt)<=ly/2);
obo   = (oby')*obx;

% Applying obstruction in optic field
g      = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(2)
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');

%% Parametrization of obstruction for rays
% np points in obstruction
no  = 5;
np  = 2^no;
th  = 2*pi/np;

TotalRays  = np;

rayH11(Nz) = OpticalRay();
rayH12(Nz) = OpticalRay();
rayH21(Nz) = OpticalRay();
rayH22(Nz) = OpticalRay();

%puntos sobre el rectangulo dada su parametrizacion, donde empezaran los
%rayos
for point_index = 1:TotalRays
  
  xj         = xt + (lx/2)*(abs(cos((point_index)*th))*cos((point_index)*th)...
             + abs(sin((point_index)*th))*sin((point_index)*th));
  yj         = yt + (ly/2)*(abs(cos((point_index)*th))*cos((point_index)*th)...
             - abs(sin((point_index)*th))*sin((point_index)*th));
  zj         = 0;

  hankeltype = 11;
  rayH11(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH11(1),point_index,hankeltype);
  hankeltype = 12;
  rayH12(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH12(1),point_index,hankeltype);
  hankeltype = 21;
  rayH21(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH21(1),point_index,hankeltype);
  hankeltype = 22;
  rayH22(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH22(1),point_index,hankeltype);
  
end


figure(3)
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
plotRays(rayH11(1),'r')

%% Physical Propagation

prop = paraxialPropagator(Kx,Kx',k,dz);
figure(4)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(Nx,length(z)); 
gy      = zeros(Nx,length(z));
W       = zeros(Nx,Nz,Nx);
%% propagation with respect to z
for z_index = 1:length(z)
% propagation with respect to z  
  %saving transversal fields
  gx(:,z_index)   = g(Nx/2+1,:);
  gy(:,z_index)   = g(:,Nx/2+1);
  % saving field for slices
  W (:,z_index,:) = g;
  %propagating field
  % it's needed correction in phase of FFT
  G = fftshift(fft2((g)));...*exp(-1j.*(u(1)).*(Kx)).*exp(-1j.*(u(1)).*(Kx'));
  %obtain new propagated field
  g = (ifft2(ifftshift(G.*prop)));...*exp(-1j.*(u(1)).*(X))*exp(-1j.*(u(1)).*(X'));
%%  Calculating propagation of Rays

  % propagation distance 
  zi   = z(z_index);
  % calculating Laguerre Parameters in zi
  HPzi = HermiteParameters(zi,InitialWaist,Wavelength,nu,mu);   
  
  % propagate all rays of H11
  HankelType = 11;
  [rayH11(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH11(z_index),...
                                          TotalRays,...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType); 
 
  HankelType = 12;                                      
  [rayH12(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH12(z_index),...
                                          TotalRays,...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType); 
                                        
  HankelType = 21;                                      
  [rayH21(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH21(z_index),...
                                          TotalRays,...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType);                                       
   
  HankelType = 22;                                       
  [rayH22(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH22(z_index),...
                                          TotalRays,...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType);                                       
%% plot propagate field                                      
  fig6 = figure(6);
  fig6.Position = [408 4 1037 973];
  plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
  title(['z = ', num2str(z(z_index)), ' of ', num2str(z(end)), ' microns'])

%% Plot propagated points of hankels
  plotRays(rayH11(z_index+1),'r')
  plotRays(rayH21(z_index+1),'y')
  plotRays(rayH12(z_index+1),'m')                                         
  plotRays(rayH22(z_index+1),'b')
  pause(0.01)
end


xslice = [0,x(Nx/2)]; 
yslice = [x(Nx/2)]; 
zslice = 1;

figure(7)
slice(z,x,y,abs(W),xslice,yslice,zslice)
shading interp
alpha(0.7); %axis off
set(gcf,'Color',[1 1 1])
colormap(mapgreen)

daspect([1.5,.1,.1]) 
axis tight 
view(-38.5,16) 
camzoom(1.4) 
camproj perspective
axis off

hold on

plotRaysPropagated(rayH11,rayH22,Nz);
