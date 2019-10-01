
%% add path for classes and functions
addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam

%-------------------- indices of Laguerre Gaussian Beams -----------------%
nu      = 6;
mu      = 0;

% Physical parameters [microns]
LaguerreInitialWaist = 1000;
InitialWaist         = LaguerreInitialWaist/sqrt(2*(2*nu+mu+1));
Wavelength           = 0.6328;

% Calculating gaussian parameters
GP                   = GaussianParameters(0,InitialWaist,Wavelength);

k                    = GP.k;
RayleighDistance     = GP.RayleighDistance;

% normalized parameters
% 
% lamb    = pi;
% wo      = 1;
% sigmaLo = wo*sqrt(2*(2*nu+mu+1)); 
% k       = 2*pi/lamb;
% zo      = k*wo^2/2;
% %------------------------ sampling of vectors ----------------------------%
%First we estimate samplig in z-direction with propagation distance 
% z-direction
Dz     = RayleighDistance;        % z-window (propagation distance)
Nz     = 2^5;                     % number of points in z-direction
dz     = Dz/Nz;                   % Resolution in z
z      = 0:dz:Dz;                 % z-vector z of propagation 

% waist of Laguerre Gauss Beam until z-propagation
MaxLaguerreWaist = LaguerreBeam.waistLaguerre(z(end),InitialWaist,RayleighDistance,nu,mu);

%Second we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
N      =  2^10;                   % Number of points in x,y axis
n      = -N/2+.05:N/2-1+.05;      % vector with N-points with resolution 1
Dx     = 2*MaxLaguerreWaist;     % Size of window 
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
LG  = LaguerreBeam(X,Y,0,InitialWaist,Wavelength,nu,mu);
% Optic Field to propagate 
g   = LG.OpticalField;
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
lo      = LaguerreInitialWaist/4.3; % size of obstruction in terms of waist of Laguerre
xt      = 0;                        % traslation of obstruction in x-axis
yt      = 0;                        % traslation of onstruction in y-axis
[~,rho] = cart2pol(X-xt,X'-yt);     % Convert this in polar coordinates
obo     = double(rho<=lo);          % Create Obstruction   
clear rho                           % Clean Matrix of polar coordinates
% Applying obstruction in optic field
g       = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(20)
pcolor(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)


%% ----------------------- Ray tracing (rx,z=0)  ----------------------- %%

pn        =  50;          % Number of rays
rayH1(pn) = OpticalRay;
rayH2(pn) = OpticalRay;

for jj = 1:pn
    % Cartersian coordinates of point in circunference of obstruction
    xi              = xt+lo*cos(jj*(2*pi)/(pn)); 
    yi              = yt+lo*sin(jj*(2*pi)/(pn));
    
    
    % saving cartesian coordinates of point in circunference in each ray
    rayH1(jj).xCoordinate(1) = xi; 
    rayH1(jj).yCoordinate(1) = yi;
    rayH2(jj).xCoordinate(1) = xi; 
    rayH2(jj).yCoordinate(1) = yi;
    
    [rayH1(jj)] = getLaguerreSlopes(rayH1(jj),x,y,z,...
                                              dx,dx,dz,...
                                              xi,yi,0,...
                                              InitialWaist,Wavelength,nu,mu,1);
    
    [rayH2(jj)] = getLaguerreSlopes(rayH1(jj),x,y,z,...
                                              dx,dx,dz,...
                                              xi,yi,0,...
                                              InitialWaist,Wavelength,nu,mu,2);
end

% Initial Field with rays in this init conditions
figure(3)
pcolor(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)
hold on
for jj=1:pn
    plot(rayH1(jj).xCoordinate(1)/(sqrt(2)*InitialWaist),rayH1(jj).yCoordinate(1)/(sqrt(2)*InitialWaist)...
        ,'.','MarkerSize',10,'LineWidth',2,'color','r')
end
hold off


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

    set(gca,'un','n','pos',[0,0,1,1])
    imagesc(x/(sqrt(2)*InitialWaist),x/(sqrt(2)*InitialWaist),abs(g).^2)
    axis off
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    drawnow 
    hold on
    % propagated points of H1 and H2
    for jj=1:pn
        plot(rayH1(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
             rayH1(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',20,'LineWidth',2,'color','r')
        plot(rayH2(jj).xCoordinate(ii-1)/(sqrt(2)*InitialWaist),...
             rayH2(jj).yCoordinate(ii-1)/(sqrt(2)*InitialWaist),'.','MarkerSize',20,'LineWidth',2,'color','y')
    end
   hold off
   pause(1)

    %------------------------ Calculating Rays ---------------------------%   
    % Given z we find z+dz and new position x+dx and y+dy with help of
    % slopes

    for jj=1:pn
        % step of ray in z-direction, equation or ray r(z) = mrz*dz+r(z-1)
        % dependece of last point, new vector
        rayH1(jj).xCoordinate(ii) = rayH1(jj).xCoordinate(ii-1) + (1/rayH1(jj).zxSlope)*dz;
        rayH1(jj).yCoordinate(ii) = rayH1(jj).yCoordinate(ii-1) + (1/rayH1(jj).zySlope)*dz;
      
        rayH2(jj).xCoordinate(ii) = rayH2(jj).xCoordinate(ii-1) + (1/rayH2(jj).zxSlope)*dz;
        rayH2(jj).yCoordinate(ii) = rayH2(jj).yCoordinate(ii-1) + (1/rayH2(jj).zySlope)*dz;
  
        % ------------------ Calculating Slopes ------------------------- %
        % point for H1
        xi = rayH1(jj).xCoordinate(ii);
        yi = rayH1(jj).yCoordinate(ii);
        zi = z(ii);

        [rayH1(jj)] = getLaguerreSlopes(rayH1(jj),x,y,z,...
                                                  dx,dx,dz,...
                                                  xi,yi,zi,...
                                                  InitialWaist,Wavelength,nu,mu,1);
        % point of H2 
        xi = rayH2(jj).xCoordinate(ii);
        yi = rayH2(jj).yCoordinate(ii);
        zi = z(ii);                                    
        [rayH2(jj)] = getLaguerreSlopes(rayH2(jj),x,y,z,...
                                                  dx,dx,dz,...
                                                  xi,yi,zi,...
                                                  InitialWaist,Wavelength,nu,mu,2);
    end
    %------------------------ End calculating rays -----------------------%   
    %propagating field
    G =fftshift(fft2(g));
    figure(7)
    imagesc(angle(G))
    %G = G.*exp(1i*pi*50);
    figure(8)
    imagesc(angle(G))
    %obtain new propagated field
    g=ifft2(fftshift(G.*prop));....*exp(-(X.^50+Y.^50)./Dx.^50);
    %saving transversal fields
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);
%     
end

