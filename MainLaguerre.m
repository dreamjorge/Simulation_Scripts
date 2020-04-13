
%% add path for classes and functions
addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam

%-------------------- indices of Laguerre Gaussian Beams -----------------%
nu = 19;
mu = 4;

% Calculating Laguerre Parameters in z = 0

% Physical parameters [microns]
InitialWaist         = 100;%179*2;
Wavelength           = 0.6328;

% normalized parameters
% InitialWaist         = 1;
% Wavelength           = pi;
% Estimating parameter of laguerre in z=0
PropagationDistance  = 0;
LP                   = LaguerreParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu);

% Initial Waist of Laguerre Beam
LaguerreInitialWaist = LP.LaguerreWaist;

% %------------------------ sampling of vectors ----------------------------%
%First we estimate sampling in z-direction with propagation distance 
% z-direction
Dz = LP.RayleighDistance;    % z-window (propagation distance)
Nz = 2^5;                    % number of points in z-direction
dz = Dz/Nz;                  % Resolution in z
z  = 0:dz:Dz;                % z-vector z of propagation 

%After sampling z vector, we estimage sampling in x,y-direction in terms of
%waist of max waist Laguerre Gauss Beam until max z-propagation

MaxLaguerreWaist = LaguerreParameters.waistFunction(z(end),InitialWaist,LP.RayleighDistance,nu,mu);

N   = 2^9;                   % Number of points in x,y axis
n   = -N/2+.05:N/2-1+.05;    % vector with N-points with resolution 1
Dx  = 2.2*MaxLaguerreWaist;  % Size of window 
dx  = Dx/N;                  % Resolution
x   = n*dx;                  % Vector with dimentions
y   = x;                     % same for y
[X] = meshgrid(x,y);         % Matrix for x,y

% transformation of coordinates
[dth,dr]                      = cart2pol(dx,dx');
[thetaCoordinate,rCoordinate] = cart2pol(x,x');
[ThetaCoordinate,RCoordinate] = cart2pol(X,X');

%diferential vector in cylindrical coordinates
difr = [dx,dx,dz];

% Last we estimate vectors of frequency for Fourier Transforms associated
% with x,y

Du   = 1/dx;                 % Size of window 
du   = 1/Dx;                 % Resolution
u    = n*du;                 % freq vector with dimentions
[U]  = meshgrid(u);          % freq Matrix
%kx,ky vectors
kx   = 2*pi*u;               % angular freq vector
[Kx] = meshgrid(kx);         % angular freq matrix

%% ----------------------- Laguerre Gauss in z = 0 --------------------- %%
PropagationDistance  = 0;
LG                   = LaguerreBeam(RCoordinate,ThetaCoordinate,PropagationDistance,InitialWaist,Wavelength,nu,mu);

% Plot of Field
figure(1)
pcolor(x, x, abs(LG.OpticalField).^2)
axis square
shading flat
colormap(mapgreen)
plotCircle(0,0,LP.LaguerreWaist);
xlabel('$x \left[ microns \right]$','Interpreter','latex') 
ylabel('$y \left[ microns \right]$','Interpreter','latex') 

% Optic Field to propagate 
g   = LG.OpticalField;
% Max Peak
pxy = max(max(g));

%% ----------------- Obstruction on Lagurre in z = 0 ------------------- %%
PropagationDistance = 0;
lo                  = (LaguerreParameters.waistFunction(...
                       PropagationDistance,InitialWaist,LP.RayleighDistance,nu,mu))/4.3;              % size of obstruction in terms of waist of Laguerre
xt                  = 0;                                                                              % traslation of obstruction in x-axis
yt                  = 0;                                                                              % traslation of onstruction in y-axis
[~,rho]             = cart2pol(X-xt,X'-yt);                                                           % Convert this in polar coordinates
obo                 = double(rho<=lo);                                                                % Create Obstruction   
clear rho                                                                                             % Clean Matrix of polar coordinates
% Applying obstruction in optic field
g                   = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(20)
pcolor(x,x,abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)


%% ----------------------- Ray tracing (rx,z=0)  ----------------------- %%

pn           =  8;          % Number of rays
rayH1(pn)    = OpticalRay;
rayH2(pn)    = OpticalRay;
SlopesH1(pn) = Slopes;
SlopesH2(pn) = Slopes;
% temp th for estimate cross in origin 
thTempPrev  = zeros(1,pn);
thTemp      = zeros(1,pn);
rObst       = zeros(1,pn);
rAcum       = zeros(1,pn);
rTempPrev   = zeros(1,pn);
rTempActual = zeros(1,pn);
rayTemp     = OpticalRay;

for point_index = 1:pn
    % Cartersian coordinates of point in circunference of obstruction
    xi = xt+lo*cos(point_index*(2*pi)/(pn)); 
    yi = yt+lo*sin(point_index*(2*pi)/(pn));
    zi = 0;
    %estimate radii of point (x,y) to origin
    [thi,ri] = cart2pol(xi,yi);
  
    %for iterative process in propagation 
    rTempPrev(point_index) = rObst(point_index);
    
    % saving cartesian coordinates of point in circunference in each ray
    rayH1(point_index).xCoordinate(1)     = xi; 
    rayH1(point_index).yCoordinate(1)     = yi;
    rayH1(point_index).zCoordinate(1)     = zi;
    rayH2(point_index).xCoordinate(1)     = xi; 
    rayH2(point_index).yCoordinate(1)     = yi;
    rayH2(point_index).zCoordinate(1)     = zi;
    rayH1(point_index).thetaCoordinate(1) = ri;
    rayH1(point_index).rCoordinate(1)     = thi;
    rayH2(point_index).thetaCoordinate(1) = ri;
    rayH2(point_index).rCoordinate(1)     = thi;    
    
    rayTemp.xCoordinate          = xi;
    rayTemp.yCoordinate          = yi;
    rayTemp.zCoordinate          = zi;
    rayTemp.rCoordinate          = ri;
    rayTemp.thetaCoordinate      = thi;
    
    HankelType  = 1;
    [SlopesH1(point_index)] = HankelLaguerre.getLaguerreSlopes(rayTemp,rCoordinate,thetaCoordinate,z,...
                                                      difr,...
                                                      LP,HankelType);
    HankelType  = 2;
    [SlopesH2(point_index)] = HankelLaguerre.getLaguerreSlopes(rayTemp,rCoordinate,thetaCoordinate,z,...
                                                      difr,...
                                                      LP,HankelType);
end

% Initial Field with rays in this init conditions
figure(3)
pcolor(x,x,abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)
hold on
for point_index=1:pn
    plot(rayH1(point_index).xCoordinate(1),rayH1(point_index).yCoordinate(1)...
        ,'.','MarkerSize',10,'LineWidth',2,'color','r')
end
hold off


%%  ----------------------- Physical Propagation ------------------------ %
% paraxial propagator 
k       = LP.k;
prop    = exp(1i*dz*(Kx.^2+(Kx').^2)/(2*k));
figure(5)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(N,length(z)); 
gy      = zeros(N,length(z));
% Save field in z = 0 
gx(:,1) = g(N/2+1,:);
gy(:,1) = g(:,N/2+1);


for z_index = 2:length(z) % propagation with respect to z
  
    % field before propagation i.e in z(ii-1)
    pxyz         = g(1,1);
    g(1,1)       = pxy;
    fig          = figure(6);

    fig.Position = [ 239 135 1354 733];
    set(gca,'un','n','pos',[0,0,1,1])
    imagesc(x,x,abs(g).^2)
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    drawnow 
    hold on
    
    % propagated points of H1 and H2 in iteration before
    for point_index=1:pn
        plot(rayH1(point_index).xCoordinate(z_index-1),...
             rayH1(point_index).yCoordinate(z_index-1),'.','MarkerSize',20,'LineWidth',2,'color','r')
        plot(rayH2(point_index).xCoordinate(z_index-1),...
             rayH2(point_index).yCoordinate(z_index-1),'.','MarkerSize',20,'LineWidth',2,'color','y')
    end
    hold off
    pause(1)

    %------------------------ Calculating Rays ---------------------------%   
    % Given z we find z+dz and new position x+dx and y+dy with help of
     % slopes
     
    %propagation distance 
    zi          = z(z_index);
    for point_index = 1:pn
        % step of ray in z-direction, equation or ray r(z) = mrz*dz+r(z-1)
        % dependece of last point, new vector
 
        rayH1(point_index).xCoordinate(z_index) = rayH1(point_index).xCoordinate(z_index-1) + (1./SlopesH1(point_index).zx)*dz;
        rayH1(point_index).yCoordinate(z_index) = rayH1(point_index).yCoordinate(z_index-1) + (1./SlopesH1(point_index).zy)*dz;

        rayH2(point_index).xCoordinate(z_index) = rayH2(point_index).xCoordinate(z_index-1) + (1./SlopesH2(point_index).zx)*dz;
        rayH2(point_index).yCoordinate(z_index) = rayH2(point_index).yCoordinate(z_index-1) + (1./SlopesH2(point_index).zy)*dz;        
        
        %point of H1 for iteration 
        xi = rayH1(point_index).xCoordinate(z_index);
        yi = rayH1(point_index).yCoordinate(z_index);
        
        [thi,ri]                = cart2pol(xi,yi);
        rayTemp.rCoordinate     = ri;
        rayTemp.thetaCoordinate = thi;        
        raytemp.zCoordinate     = zi;
        
        HankelType  = 1;
        [SlopesH1(point_index)] = HankelLaguerre.getLaguerreSlopes(rayTemp,...
                                                                   rCoordinate,thetaCoordinate,z,...
                                                                   difr,...
                                                                   LP,HankelType); 
        %point of H2 for iteration
        xi = rayH2(point_index).xCoordinate(z_index);
        yi = rayH2(point_index).yCoordinate(z_index);
        
        [thi,ri]                = cart2pol(xi,yi);
        rayTemp.rCoordinate     = ri;
        rayTemp.thetaCoordinate = thi;        
        raytemp.zCoordinate     = zi;
                
        %----------condition for cross around origin
        rTempActual(point_index) = ri;
        drTemp                   = abs( rTempPrev(point_index) - rTempActual(point_index) );
        rAcum(point_index)       = (rAcum(point_index) + drTemp);
        
        if ( rAcum(point_index)<rObst(point_index) )
          
          HankelType  = 2;
          [SlopesH2(point_index)] = HankelLaguerre.getLaguerreSlopes(rayTemp,...
                                                                     rCoordinate,thetaCoordinate,z,...
                                                                     difr,...
                                                                     LP,HankelType);
        else % change to H1
           
          HankelType  = 1;
          [SlopesH2(point_index)] = HankelLaguerre.getLaguerreSlopes(rayTemp,...
                                                                     rCoordinate,thetaCoordinate,z,...
                                                                     difr,...
                                                                     LP,HankelType);                                    
        end
        
        rTempPrev(point_index) = rTempActual(point_index);
                                                     
    end
    %------------------------ End calculating rays -----------------------%   
    %propagating field
    % it's needed correction in phase of FFT
    G = fftshift(fft2(ifftshift(g)));...*exp(-1j.*(u(1)).*(Kx)).*exp(-1j.*(u(1)).*(Kx'));
    %obtain new propagated field
    g = fftshift(ifft2(ifftshift(G.*prop)));...*exp(-1j.*(u(1)).*(X))*exp(-1j.*(u(1)).*(X'));
    figure(7)
    imagesc(angle(G))
    title('Angle of Optical Field')
    %G = G.*exp(1i*pi*50);
    figure(8)
    imagesc(angle(g))
    title('Angle of Optical Field FT')
    %obtain new propagated field
    %saving transversal fields
    gx(:,z_index) = g(N/2+1,:);
    gy(:,z_index) = g(:,N/2+1);
    
    figure(9)
    H1 = XLaguerreBeam(RCoordinate,ThetaCoordinate,zi,InitialWaist,Wavelength,nu,mu);
    imagesc(angle(H1.OpticalField))
    title('Angle of Exact H1 Optical Field')
    
    figure(10)
    H2 = HankelLaguerre(RCoordinate,ThetaCoordinate,zi,InitialWaist,Wavelength,nu,mu,2);
    imagesc(angle(H2.OpticalField))
    title('Angle of Exact H2 Optical Field')
%     
%     
end

