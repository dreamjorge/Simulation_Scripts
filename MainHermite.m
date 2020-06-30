%%                      Main Hermite
% Ugalde-Ontiveros J.A. 
%% add path for classes and functions
clear all
addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam
GenerateVideo = 'NO';
%% indices of Hermite Gaussian Beams 
nu = 12;
mu = 11;

%% Physical parameters [microns]
InitialWaist          = 100;
Wavelength            = 0.6328;

HPz0                  = HermiteParameters(0,InitialWaist,Wavelength,nu,mu);

k                     = HPz0.k;
RayleighDistance      = HPz0.RayleighDistance;

HermiteInitialWaistX  = HPz0.HermiteWaistX;
HermiteInitialWaistY  = HPz0.HermiteWaistY;
HermiteInitialWaist   = HPz0.HermiteWaist;

%% Normalized parameters
% 
% Wavelength           = pi;
% InitialWaist         = 1;
% HermiteInitialWaist  = InitialWaist*sqrt(nu+mu+1);
% HPz0                 = HermiteParameters(0,InitialWaist,Wavelength,nu,mu);
% k                    = HPz0.k;
% RayleighDistance     = HPz0.RayleighDistance;

%% sampling of vectors 
%First, we estimate samplig in z-direction with propagation distance 
% z-direction
Dz    = 1.5*RayleighDistance;   % z-window (propagation distance)
Nz    = 2^7;                  % number of points in z-direction
dz    = Dz/Nz;                % Resolution in z
nz    = 0:Nz-1;               % vector with N-points with resolution 1
z     = nz*dz;                % z-vector z of propagation 

% waist of Hermite Gauss Beam until z-propagation
MaxHermiteWaist = HermiteParameters.getWaist(z(end),InitialWaist,RayleighDistance,nu,mu);

%Second, we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
Nx    =  2^9;                % Number of points in x,y axis
n     = -Nx/2:Nx/2-1;         % vector with N-points with resolution 1
Dx    = MaxHermiteWaist;      % Size of window 
dx    = Dx/Nx;                % Resolution
x     = n*dx;                 % Vector
y     = x;
[X,Y] = meshgrid(x,y);

%Last, we estimate vectors of frequency for Fourier Transforms
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

% Changing parameters to vector z
HPz.zCoordinate = z;  

% Optic Field to propagate 
g     = HGB.OpticalFieldHermite;
% Plot of Function
figure(1)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),abs(g),mapgreen,'');
saveas(gcf,'HermiteBeam.png')


%% Hankels 
H11 = HankelHermite(X,Y,HPz0,11);
H12 = HankelHermite(X,Y,HPz0,12);
H21 = HankelHermite(X,Y,HPz0,21);
H22 = HankelHermite(X,Y,HPz0,22);
figure(11)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(H11.OpticalField),parula,'');
title('Phase of Hermite Hankel [1,1]')
colorbar
saveas(gcf,'HH11.png')
figure(12)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(H12.OpticalField),parula,'');
title('Phase of Hermite Hankel [1,2]')
colorbar
saveas(gcf,'HH12.png')
figure(13)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(H21.OpticalField),parula,'');
title('Phase of Hermite Hankel [2,1]')
colorbar
saveas(gcf,'HH21.png')
figure(14)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(H22.OpticalField),parula,'');
title('Phase of Hermite Hankel [2,2]')
colorbar
saveas(gcf,'HH22.png')

figure(15)
sum1122 = (H11.OpticalField + H22.OpticalField);
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(sum1122),parula,'');
title('Phase of (Hermite Hankel [1,1]+[2,2]')
colorbar
saveas(gcf,'angleHH11plus22.png')
sum112212 = sum1122 + H12.OpticalField;
figure(16)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(sum112212),parula,'');
title('Phase of (Hermite Hankel [1,1]+[2,2]+[1,2]')
colorbar
saveas(gcf,'angleHH11plus22plus12.png')
sum1112 = (H11.OpticalField + H12.OpticalField);
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),angle(sum1112),parula,'');
title('Phase of (Hermite Hankel [1,1]+[2,2]')
colorbar
saveas(gcf,'angleHH11plus12.png')
%%

      WaistGauss = HPz0.Waist;
      PhiPhase   = HPz0.PhiPhase;

      [Hx,NHx] = ...
      HermiteBeam.hermiteSolutions(nu...
                                  ,(sqrt(2)./WaistGauss).*x);
      
                                
      GaussX   = GaussianBeam(x,HPz0).OpticalField;
      
      %Hx       = Hx .*GaussX.*exp(1i*PhiPhase);
      %NHx      = NHx.*GaussX.*exp(1i*PhiPhase);

      H1x      = Hx+1i*NHx;
      
      figure(16)
      plot(x/(HPz0.Waist),abs(Hx),'LineWidth',1.6)
      hold on
      plot(x/(HPz0.Waist),abs(NHx),'--','LineWidth',1.6)
      plot(x/(HPz0.Waist),abs(H1x),'-.','LineWidth',1.6)
      hold off
      xlim([-HPz0.HermiteWaistX/(2.3*HPz0.Waist),HPz0.HermiteWaistX/(2.3*HPz0.Waist)])
      ylim([0 60])
      xlabel('$\xi$','Interpreter','latex')
      legend('Hermite','NHermite',' Hankel Hermite')
      title('Solutions of Hermite Equation of integer number equal 16')
      saveas(gcf,'hankelHermiteSolutions.png')


%% Obstruction on Hermite in z = 0

% Obstruction
lx    = HPz0.HermiteWaist/3.8;
ly    = HPz0.HermiteWaist/4.5;
xt    = 0;...HPz0.HermiteWaist/5;
yt    = xt;
obx   = double(abs(x-xt)<=lx/2);
oby   = double(abs(x-yt)<=ly/2);
obo   = (oby')*obx;

% Applying obstruction in optic field
g      = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(2)
plotOpticalField(x/(HPz0.Waist),x/(HPz0.Waist),abs(g),mapgreen,'');
saveas(gcf,'HermiteBeamWithObstruction.png')
%% Parametrization of obstruction for rays
% Total points/rays in obstruction
no        = 5;
TotalRays = 2^no;
th        = 2*pi/TotalRays;

% Optical Ray of size Nz, i.e size of direction of z-coordinate
rayH11(Nz) = OpticalRay();
rayH12(Nz) = OpticalRay();
rayH21(Nz) = OpticalRay();
rayH22(Nz) = OpticalRay();

% it given parametrization of obstruction, we give initial conditions of
% rays
for ray_index = 1:TotalRays
  
  xj         = xt + (lx/2)*(abs(cos((ray_index)*th))*cos((ray_index)*th)...
             + abs(sin((ray_index)*th))*sin((ray_index)*th));
  yj         = yt + (ly/2)*(abs(cos((ray_index)*th))*cos((ray_index)*th)...
             - abs(sin((ray_index)*th))*sin((ray_index)*th));
  zj         = 0;

  %assign coordinates to each Hankel ray
  hankeltype = 11;
  rayH11(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH11(1),ray_index,hankeltype);
  hankeltype = 12;
  rayH12(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH12(1),ray_index,hankeltype);
  hankeltype = 21;
  rayH21(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH21(1),ray_index,hankeltype);
  hankeltype = 22;
  rayH22(1)  = assignCoordinates2CartesianRay(xj,yj,zj,rayH22(1),ray_index,hankeltype);
  
end

figure(3)
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
plotRays(rayH11(1),'r')
saveas(gcf,'HermiteBeamWithObstructionRays.png')
%% Physical Propagation

prop = paraxialPropagator(Kx,Kx',k,dz);
figure(4)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(Nx,length(z)); 
gy      = zeros(Nx,length(z));
W       = zeros(Nx,Nz,Nx);

% Videos Options for generate video
if strcmp(GenerateVideo,'YES') % this was defined 8-line
  vidObj1 = VideoWriter('HGo.avi');
  vidObj1.Quality   = 100;
  vidObj1.FrameRate = 30;
  open(vidObj1);
end
%% propagation with respect to z
for z_index = 1:length(z)
  %% plot propagate field                                      
  fig6 = figure(6);
  fig6.Position = [408 4 1037 973];
  plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
  hold on
%% Plot propagated points of hankels
  plotRays(rayH11(z_index),'r')
  plotRays(rayH21(z_index),'y')
  plotRays(rayH12(z_index),'m')                                         
  plotRays(rayH22(z_index),'c')
  title(['z = ', num2str(z(z_index)), ' of ', num2str(z(end)), ' microns'])
  drawnow
% Write video
  if strcmp(GenerateVideo,'YES')
    writeVideo(vidObj1, getframe(gca));
  end                                    
        
% propagation with respect to z  
  % saving transversal fields
  gx(:,z_index)   = g(Nx/2+1,:);
  gy(:,z_index)   = g(:,Nx/2+1);
  % saving field for slices
  W (:,z_index,:) = g;
  % propagating field
  G = fftshift(fft2((g)));
  % obtain new propagated field
  g = (ifft2(ifftshift(G.*prop)));
%%  Calculating propagation of Rays
  % propagation distance 
  zi   = z(z_index);
  % calculating Laguerre Parameters in zi
  HPzi = HermiteParameters(zi,InitialWaist,Wavelength,nu,mu);   
  
  % propagate all rays of H11
  HankelType = 11;
  [rayH11(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH11(z_index),...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType); 
 
  HankelType = 12;                                      
  [rayH12(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH12(z_index),...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType); 
                                        
  HankelType = 21;                                      
  [rayH21(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH21(z_index),...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType);                                       
   
  HankelType = 22;                                       
  [rayH22(z_index+1)] = ...
  HankelHermite.getPropagateCartesianRays(rayH22(z_index),...
                                          x,y,...
                                          dr,...
                                          HPzi,...
                                          HPz,...
                                          HankelType);                                       
                                
    pause(0.01)
end

%% plot propagate field at z(end)                                    
fig6 = figure(6);
fig6.Position = [408 4 1037 973];
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
title(['z = ', num2str(z(z_index)), ' of ', num2str(z(end)), ' microns'])

%% Plot propagated points of hankels at z(end)
plotRays(rayH11(z_index+1),'r')
plotRays(rayH21(z_index+1),'y')
plotRays(rayH12(z_index+1),'m')
plotRays(rayH22(z_index+1),'c')

if strcmp(GenerateVideo,'YES')
  writeVideo(vidObj1, getframe(gca));
  close(vidObj1);
end
%% Slices 
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

plotPropagatedRays(rayH11,rayH22);
%%
zindex = 128;
gg = W(:,zindex,:);
gg = reshape(gg,[Nx,Nx]);
figure(8)
plotOpticalField(x,x,abs(g).^1.5,mapgreen,'microns');
plotRays(rayH11(zindex),'r')
plotRays(rayH21(zindex),'y')
plotRays(rayH12(zindex),'m')                                         
plotRays(rayH22(zindex),'c')