
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
z      = 0:dz:2*Dz;                   % z-vector z of propagation 

% waist of Laguerre Gauss Beam until z-propagation
maxIndex        = max(nu,mu);
MaxHermiteWaist = HermiteBeam.waistHermite(z(end),InitialWaist,RayleighDistance,maxIndex);

%Second we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
N      =  2^8;                   % Number of points in x,y axis
n      = -N/2+.05:N/2-1+.05;     % vector with N-points with resolution 1
Dx     = 2*MaxHermiteWaist;      % Size of window 
dx     = Dx/N;                   % Resolution
x      = n*dx;                   % Vector
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

% vector r
r    = struct();
r.x  = x;
r.y  = y;
r.z  = z;

% differential of vector
dr   = struct();
dr.x = dx;
dr.y = dx;
dr.z = dz;
  
%% ----------------------- Hermite Gauss in z = 0 --------------------- %%
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

%% ----------------- Obstruction on Hermite in z = 0 ------------------- %%

%Estimating parameters of Hermite in z=0
HParameters = HermiteParameters(0,InitialWaist,Wavelength,nu,mu);

% Obstruction
lx     = HParameters.HermiteWaist/3;
ly     = lx;
xt     = 0;
yt     = 0;
obx    = double(abs(x)<=lx/2);
oby    = double(abs(x)<=ly/2);
obo    = (oby')*obx;

% Applying obstruction in optic field
g      = g.*(1-obo);
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


% parametrization of obstruction for rays
% np points in obstruction
no  = 5;
np  = 2^no;
th  = 2*pi/np;

TotalNumberRays         = np;
rayH11(TotalNumberRays) = OpticalRay;
rayH12(TotalNumberRays) = OpticalRay;
rayH21(TotalNumberRays) = OpticalRay;
rayH22(TotalNumberRays) = OpticalRay;

%puntos sobre el rectangulo dada su parametrizacion, donde empezaran los
%rayos
for jj = 1:TotalNumberRays
  
  xj = xt+(lx/2)*(abs(cos((jj)*th))*cos((jj)*th)+abs(sin((jj)*th))*sin((jj)*th));
  yj = yt+(ly/2)*(abs(cos((jj)*th))*cos((jj)*th)-abs(sin((jj)*th))*sin((jj)*th));
  
  rayH11(jj).xCoordinate(1) = xj;
  rayH11(jj).yCoordinate(1) = yj;
  
  rayH12(jj).xCoordinate(1) = xj;
  rayH12(jj).yCoordinate(1) = yj;
  
  rayH21(jj).xCoordinate(1) = xj;
  rayH21(jj).yCoordinate(1) = yj;
  
  rayH22(jj).xCoordinate(1) = xj;
  rayH22(jj).yCoordinate(1) = yj;
    
  % vector of propagation
  rp   = struct();
  rp.x = xj;
  rp.y = yj;
  rp.z = z(1); 
  
  % Estimate slopes of in rp
  HankelType   = [1,1];
  [rayH11(jj)] = HankelHermite.getHermiteSlopes(rayH11(jj),r,dr,rp,HParameters,HankelType);

  HankelType   = [1,2];
  [rayH12(jj)] = HankelHermite.getHermiteSlopes(rayH12(jj),r,dr,rp,HParameters,HankelType);

  HankelType   = [2,1];
  [rayH21(jj)] = HankelHermite.getHermiteSlopes(rayH21(jj),r,dr,rp,HParameters,HankelType);  
  
  HankelType   = [2,2];
  [rayH22(jj)] = HankelHermite.getHermiteSlopes(rayH22(jj),r,dr,rp,HParameters,HankelType);  
                               
        
end


figure(3)
pcolor(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;

hold on
for jj=1:TotalNumberRays
  plot(rayH11(jj).xCoordinate(1)/(sqrt(2)*InitialWaist),rayH11(jj).yCoordinate(1)/(sqrt(2)*InitialWaist),'.','color','r')
end
hold off
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

% set(gca,'un','n','pos',[0,0,1,1])
  imagesc(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
  colormap(mapgreen)
  set(gca,'YDir','normal')
  axis square
  title(['z = ', num2str(z(ii))])
  drawnow 
  hold on
  for jj = 1:TotalNumberRays
  
    plot(rayH11(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
         rayH11(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',10,'LineWidth',2,'color','r')

    plot(rayH12(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
         rayH12(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',10,'LineWidth',2,'color','y')

    plot(rayH21(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
         rayH21(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',10,'LineWidth',2,'color','b')

    plot(rayH22(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
         rayH22(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',10,'LineWidth',2,'color','g')      

    %% new point of propagation   
    rayH11(jj).xCoordinate(ii) = rayH11(jj).xCoordinate(ii-1) + (1/rayH11(jj).zxSlope)*dz;
    rayH11(jj).yCoordinate(ii) = rayH11(jj).yCoordinate(ii-1) + (1/rayH11(jj).zySlope)*dz;

    rayH12(jj).xCoordinate(ii) = rayH12(jj).xCoordinate(ii-1) + (1/rayH12(jj).zxSlope)*dz;
    rayH12(jj).yCoordinate(ii) = rayH12(jj).yCoordinate(ii-1) + (1/rayH12(jj).zySlope)*dz;

    rayH21(jj).xCoordinate(ii) = rayH21(jj).xCoordinate(ii-1) + (1/rayH21(jj).zxSlope)*dz;
    rayH21(jj).yCoordinate(ii) = rayH21(jj).yCoordinate(ii-1) + (1/rayH21(jj).zySlope)*dz;

    rayH22(jj).xCoordinate(ii) = rayH22(jj).xCoordinate(ii-1) + (1/rayH22(jj).zxSlope)*dz;
    rayH22(jj).yCoordinate(ii) = rayH22(jj).yCoordinate(ii-1) + (1/rayH22(jj).zySlope)*dz;
    
    %% new slopes in new point
    
    % ray for H11
    rp.x = rayH11(jj).xCoordinate(ii);
    rp.y = rayH11(jj).yCoordinate(ii);
    rp.z = z(ii); 

    % Estimate slopes of in rp
    HankelType   = [1,1];
    [rayH11(jj)] = HankelHermite.getHermiteSlopes(rayH11(jj),r,dr,rp,HParameters,HankelType);

    % ray for H12
    rp.x = rayH12(jj).xCoordinate(ii);
    rp.y = rayH12(jj).yCoordinate(ii);
    rp.z = z(ii); 
    
    % Estimate slopes of in rp
    HankelType   = [1,2];
    [rayH12(jj)] = HankelHermite.getHermiteSlopes(rayH12(jj),r,dr,rp,HParameters,HankelType);

    % ray for H21
    rp.x = rayH21(jj).xCoordinate(ii);
    rp.y = rayH21(jj).yCoordinate(ii);
    rp.z = z(ii);     
    
    % Estimate slopes of in rp
    HankelType   = [2,1];
    [rayH21(jj)] = HankelHermite.getHermiteSlopes(rayH21(jj),r,dr,rp,HParameters,HankelType);  

    % ray for H22
    rp.x = rayH22(jj).xCoordinate(ii);
    rp.y = rayH22(jj).yCoordinate(ii);
    rp.z = z(ii);    
    
    % Estimate slopes of in rp    
    HankelType   = [2,2];
    [rayH22(jj)] = HankelHermite.getHermiteSlopes(rayH22(jj),r,dr,rp,HParameters,HankelType);                                              

  end

  hold off
  pause(.01)

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
  imagesc(unwrap(angle(G)))
  %obtain new propagated field
  %saving transversal fields
  gx(:,ii)=g(N/2+1,:);
  gy(:,ii)=g(:,N/2+1);

    
%     
end

