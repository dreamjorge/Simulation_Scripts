%%                      Main Hermite
% Ugalde-Ontiveros J.A. 
%% add path for classes and functions
clear all
addpath ParaxialBeams
addpath ParaxialBeams\Addons
addpath ParaxialBeams\Addons\export_fig-master
addpath ParaxialBeams\Addons\panel-2.14
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam
GenerateVideo = 'NO';
%% indices of Hermite Gaussian Beams 
nu = 12;
mu = 11;

%% Physical parameters [microns]
InitialWaist          = 100;
Wavelength            = 0.6328;

%% Build parameters of Hermite in z=0
HermiteParametersz0  = HermiteParameters(0,InitialWaist,Wavelength,nu,mu);
%obtain parameters of beam
k                    = HermiteParametersz0.k;
RayleighDistance     = HermiteParametersz0.RayleighDistance;

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
Dz    = 0.5*RayleighDistance;     % z-window (propagation distance)
Nz    = 2^7;                      % number of points in z-direction
dz    = Dz/Nz;                    % Resolution in z
nz    = 0:Nz-1;                   % vector with N-points with resolution 1
z     = nz*dz;                    % z-vector z of propagation 

% waist of Hermite Gauss Beam until z-propagation
MaxHermiteWaist = HermiteParameters.getWaist(z(end),InitialWaist,RayleighDistance,nu,mu);

%Second, we estimate sampling in x,y-direction in terms of waist of Guassian
%Hermite Beam

% y,x-direction
Nx    =  2^10;                 % Number of points in x,y axis
n     = -Nx/2:Nx/2-1;         % vector with N-points with resolution 1
Dx    = 1.1*MaxHermiteWaist;  % Size of window 
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

%% Hankel 1D

HermiteInitialWaistX  = HermiteParametersz0.HermiteWaistX;
HermiteInitialWaistY  = HermiteParametersz0.HermiteWaistY;
HermiteInitialWaist   = HermiteParametersz0.HermiteWaist;
PhiPhase              = HermiteParametersz0.PhiPhase;
% for x-direction
InitialWaist          = HermiteParametersz0.InitialWaist;
[Hx,NHx]              = ...
HermiteBeam.hermiteSolutions(nu,(sqrt(2)./InitialWaist).*x);

GaussX                = GaussianBeam(x,HermiteParametersz0).OpticalField;
Hx                    = Hx.*GaussX.*exp(1i*PhiPhase);
NHx                   = NHx.*GaussX.*exp(1i*PhiPhase);

%trunk with gaussian
SupGaussX             = exp(-(sqrt(2)*x./(HermiteInitialWaistX)).^(50)); 
NHx                   = NHx.*SupGaussX;

H1x                   = Hx+1j*NHx;
H2x                   = Hx-1j*NHx;
%% plots for x
close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hx),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGauss1SolutionX','-png','-transparent')

close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(NHx),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$NH_n$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGauss2SolutionX','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(H1x),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$HH_n^{1}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussHankel1SolutionX','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(H2x),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$HH_n^{2}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussHanke2SolutionX','-png','-transparent')


close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hx),'LineWidth',1.6)
hold on
plot(x/(HermiteParametersz0.Waist),abs(NHx),'--','LineWidth',1.6)
plot(x/(HermiteParametersz0.Waist),abs(H1x),'-.','LineWidth',1.6)
hold off
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$','$NH_n$','$HH_n^{1}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussSolutions1X','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hx),'LineWidth',1.6)
hold on
plot(x/(HermiteParametersz0.Waist),abs(NHx),'--','LineWidth',1.6)
plot(x/(HermiteParametersz0.Waist),abs(H2x),'-.','LineWidth',1.6)
hold off
xlim([-7 7])
ylim([0 .019])
xlabel('$\xi$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$','$NH_n$','$HH_n^{2}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussSolutions2X','-png','-transparent')

%% 
% for y-direction
[Hy,NHy]  = ...
HermiteBeam.hermiteSolutions(mu,(sqrt(2)./InitialWaist).*x);

% Hermite Gauss 1D and Hankel Hermite Gauss 
PhiPhase  = HermiteParametersz0.PhiPhase;
GaussX    = GaussianBeam(x,HermiteParametersz0).OpticalField;
Hy        = Hy .*GaussX.*exp(1i*PhiPhase);
NHy       = NHy.*GaussX.*exp(1i*PhiPhase);

SupGaussY = exp(-(sqrt(2)*y./(HermiteInitialWaistY)).^50);
NHy       = NHy.*SupGaussY;
H1y       = Hy+1i*NHy;
H2y       = Hy-1i*NHy;

%% plots for y
close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hy),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGauss1SolutionY','-png','-transparent')

close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(NHy),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$NH_n$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGauss2SolutionY','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(H1y),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$HH_n^{1}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussHankel1SolutionY','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(H2y),'LineWidth',1.6)
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$HH_n^{2}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussHanke2SolutionY','-png','-transparent')


close(figure(1))
fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hy),'LineWidth',1.6)
hold on
plot(x/(HermiteParametersz0.Waist),abs(NHy),'--','LineWidth',1.6)
plot(x/(HermiteParametersz0.Waist),abs(H1y),'-.','LineWidth',1.6)
hold off
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$','$NH_n$','$HH_n^{1}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussSolutions1Y','-png','-transparent')

fig1          = figure(1);
fig1.Position = [609 206 945 348];
plot(x/(HermiteParametersz0.Waist),abs(Hy),'LineWidth',1.6)
hold on
plot(x/(HermiteParametersz0.Waist),abs(NHy),'--','LineWidth',1.6)
plot(x/(HermiteParametersz0.Waist),abs(H2y),'-.','LineWidth',1.6)
hold off
xlim([-7 7])
ylim([0 .019])
xlabel('$\eta$','Interpreter','latex')
ylabel('Amplitude')
leg1=legend('$H_n$','$NH_n$','$HH_n^{2}$');
set(leg1,'Interpreter','latex');
export_fig('ImagesThesis/hankelHermiteGaussSolutions2Y','-png','-transparent')



%% Hermite Gauss in z = 0
HGB   = HermiteBeam(X,Y,HermiteParametersz0);
% copy of parameters
HPz   = copy(HermiteParametersz0);
% Changing parameters to vector z
HPz.zCoordinate = z;  
% Optic Field to propagate 
g     = HGB.OpticalFieldHermite;
% Plot of Function
close(figure(3))
figure(3)

plotOpticalField(x/(HermiteParametersz0.Waist),x/(HermiteParametersz0.Waist),abs(g),mapgreen,'');
export_fig('ImagesThesis/HermiteGaussz0','-png','-transparent')


%% Hankels 2D
H11 = (H1y')*(H1x);
H12 = (H2y')*(H1x);
H21 = (H1y')*(H2x);
H22 = (H2y')*(H2x);

% g = H11.*H22;
%% phases
close(figure(2))
fig2 = figure(2);
fig2.Position = [680 273 721 705];
ha = tight_subplot(2,2,[.01 .03],[.06 .05],[.1 .09]);
axes(ha(1))
plotOpticalField(x/InitialWaist,x/InitialWaist,(angle(H22)),parula,'');
title('$\psi_{n,m}^{1,1}$','interpreter','latex','FontSize',10)
ha(1).XAxis.Visible = 'off';
axes(ha(2))
plotOpticalField(x/InitialWaist,x/InitialWaist,(angle(H21)),parula,'');
title('$\psi_{n,m}^{1,2}$','interpreter','latex','FontSize',10)
ha(2).XAxis.Visible = 'off';
ha(2).YAxis.Visible = 'off';
%ha(2).YAxisLocation = 'right';
axes(ha(3))
plotOpticalField(x/InitialWaist,x/InitialWaist,(angle(H12)),parula,'');
title('$\psi_{n,m}^{2,1}$','interpreter','latex','FontSize',10)
%ha(3).YAxis.Visible = 'off';
axes(ha(4))
plotOpticalField(x/InitialWaist,x/InitialWaist,(angle(H11)),parula,'');
ha(4).YAxisLocation = 'right';
ha(4).YAxis.Visible = 'off';
title('$\psi_{n,m}^{2,2}$','interpreter','latex','FontSize',10)
sgtitle('Phase of Hankels')
cb = colorbar;
cb.Location = 'eastoutside';
cb.Position =  [0.9490    0.0598    0.0255    0.8755];
export_fig('ImagesThesis/PhasesHankels','-png','-transparent')

%% abs of hankels

close(figure(3))
fig3 = figure(3);
fig3.Position = [680 164 721 814];
ha = tight_subplot(2,2,[.01 .01],[.05 .01],[.1 .01]);
axes(ha(1))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21),mapgreen,'');
ha(1).XAxis.Visible = 'off';
ha(1).YAxisLocation = 'left';
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}|$','interpreter','latex','FontSize',18)
set(gca,'FontSize',18);
axes(ha(2))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H11+H22),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
ha(2).XAxis.Visible = 'off';
ha(2).YAxis.Visible = 'off';
ha(2).YAxisLocation = 'right';
set(gca,'FontSize',18);
axes(ha(3))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21+H11),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
set(gca,'FontSize',18);
ha(3).YAxisLocation = 'left';
axes(ha(4))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21+H11+H12),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}+\psi_{n,m}^{1,2}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
ha(4).YAxisLocation = 'right';
ha(4).YAxis.Visible = 'off';
set(gca,'FontSize',18);
% sgtitle('Superposition of Hankels','FontSize',18)
%saveas(gcf,'SuperpositionOfHankels.png')
export_fig('ImagesThesis/SuperpositionOfHankels','-png','-transparent')
%%
close(figure(4))
fig4 = figure(4);
fig4.Position = [680 164 721 814];
ha = tight_subplot(2,2,[.01 .01],[.05 .01],[.1 .01]);
axes(ha(1))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H21+H11),mapgreen,'');
ha(1).XAxis.Visible = 'off';
ha(1).YAxisLocation = 'left';
title('$|\psi_{n,m}^{2,1}+\psi_{n,m}^{1,1}|$','interpreter','latex','FontSize',18)
set(gca,'FontSize',18);
axes(ha(2))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H21+H12),mapgreen,'');
title('$|\psi_{n,m}^{2,1}+\psi_{n,m}^{1,2}|$','interpreter','latex','FontSize',18)
ha(2).XAxis.Visible = 'off';
ha(2).YAxis.Visible = 'off';
ha(2).YAxisLocation = 'right';
set(gca,'FontSize',18);
axes(ha(3))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H12+H21),mapgreen,'');
title('$|\psi_{n,m}^{1,2}+\psi_{n,m}^{2,1}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
set(gca,'FontSize',18);
ha(3).YAxisLocation = 'left';
axes(ha(4))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21+H11+H12),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}+\psi_{n,m}^{1,2}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
ha(4).YAxisLocation = 'right';
ha(4).YAxis.Visible = 'off';
set(gca,'FontSize',18);
% sgtitle('Superposition of Hankels','FontSize',18)
%saveas(gcf,'SuperpositionOfHankels.png')
export_fig('ImagesThesis/SuperpositionOfHankels2','-png','-transparent')

%% abs of combinations of hankels

close(figure(4))
fig4 = figure(4);
fig4.Position = [537 535 1067 443];
ha = tight_subplot(1,3,[.01 .01],[.05 .1],[.1 .07]);

axes(ha(1))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H11),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
ha(2).YAxisLocation = 'left';
set(gca,'FontSize',18);
axes(ha(2))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21+H11),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
set(gca,'FontSize',18);
ha(2).YAxis.Visible = 'off';
axes(ha(3))
plotOpticalField(x/InitialWaist,x/InitialWaist,abs(H22+H21+H11+H12),mapgreen,'');
title('$|\psi_{n,m}^{1,1}+\psi_{n,m}^{2,1}+\psi_{n,m}^{1,2}+\psi_{n,m}^{2,2}|$','interpreter','latex','FontSize',18)
ha(3).YAxisLocation = 'right';
set(gca,'FontSize',18);
sgtitle('Superposition of Hankels','FontSize',18)
%saveas(gcf,'SuperpositionOfHankels.png')
export_fig('SuperpositionOfHankels2','-png','-transparent')


%% Obstruction on Hermite in z = 0

% Obstruction
lx    = HermiteParametersz0.HermiteWaist/3.8;
ly    = HermiteParametersz0.HermiteWaist/4.5;
xt    = 0;...HPz0.HermiteWaist/5;
yt    = xt;
obx   = double(abs(x-xt)<=lx/2);
oby   = double(abs(x-yt)<=ly/2);
obo   = (oby')*obx;

g11o  = H11.*(1-obo);
g12o  = H12.*(1-obo);
g21o  = H21.*(1-obo);
g22o  = H22.*(1-obo);
g11   = H11;
g22   = H22;
g21   = H21;

% Applying obstruction in optic field
go    = g.*(1-obo);
%Ploting Laguerre with obstruction
close(figure(2))
figure(2)
plotOpticalField(x/(HermiteParametersz0.Waist),x/(HermiteParametersz0.Waist),abs(go),mapgreen,'');
export_fig('ImagesThesis/HermiteBeamWithObstruction','-png','-transparent')
%% Parametrization of obstruction for rays
% Total points/rays in obstruction
no        = 4;
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
close(figure(3))
figure(3)
plotOpticalField(x,x,abs(go),mapgreen,'microns');
plotRays(rayH11(1),'r',1)
export_fig('ImagesThesis/HermiteBeamWithObstructionRays','-png','-transparent')
%% Physical Propagation

prop = paraxialPropagator(Kx,Kx',k,dz);
figure(4)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(Nx,length(z)); 
gy      = zeros(Nx,length(z));
W       = zeros(Nx,Nz,Nx);
W11o    = zeros(Nx,Nz,Nx);
W22o    = zeros(Nx,Nz,Nx);
W12o    = zeros(Nx,Nz,Nx);
W21o    = zeros(Nx,Nz,Nx);
W11     = zeros(Nx,Nz,Nx);
W22     = zeros(Nx,Nz,Nx);
Wo      = zeros(Nx,Nz,Nx);
W21     = zeros(Nx,Nz,Nx);
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
  plotRays(rayH11(z_index),'r',1)
  plotRays(rayH21(z_index),'y',1)
  plotRays(rayH12(z_index),'m',1)                                         
  plotRays(rayH22(z_index),'c',1)
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
  W    (:,z_index,:) = g;
  W11o (:,z_index,:) = g11o;
  W22o (:,z_index,:) = g22o;
  W12o (:,z_index,:) = g12o;
  W21o (:,z_index,:) = g21o;
  Wo   (:,z_index,:) = go;
  W11  (:,z_index,:) = g11;
  W22  (:,z_index,:) = g22;
  W21  (:,z_index,:) = g21;
  % propagating field
  G    = fftshift(fft2((g)));
  Go   = fftshift(fft2((go)));
  G11  = fftshift(fft2((g11)));
  G22  = fftshift(fft2((g22)));
  G21  = fftshift(fft2((g21)));
  G11o = fftshift(fft2((g11o)));
  G22o = fftshift(fft2((g22o)));
  G12o = fftshift(fft2((g12o)));
  G21o = fftshift(fft2((g21o)));
  
  % obtain new propagated field
  go   = (ifft2(ifftshift(Go.*prop)));
  g    = (ifft2(ifftshift(G.*prop)));
  g22  = (ifft2(ifftshift(G22.*prop)));
  g21  = (ifft2(ifftshift(G21.*prop)));
  g11  = (ifft2(ifftshift(G11.*prop)));
  g11o = (ifft2(ifftshift(G11o.*prop)));
  g22o = (ifft2(ifftshift(G22o.*prop)));
  g12o = (ifft2(ifftshift(G12o.*prop)));
  g21o = (ifft2(ifftshift(G21o.*prop)));
  
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
plotRays(rayH11(z_index+1),'r',1)
plotRays(rayH21(z_index+1),'y',1)
plotRays(rayH12(z_index+1),'m',1)
plotRays(rayH22(z_index+1),'c',1)

if strcmp(GenerateVideo,'YES')
  writeVideo(vidObj1, getframe(gca));
  close(vidObj1);
end

%% transversal field
close(figure(1))
figure(1);
imagesc(z/RayleighDistance,x/InitialWaist,abs(gy).^2)
xlabel('$z$','Interpreter','latex')
ylabel('$\eta$','Interpreter','latex')
colormap(mapgreen)


%% Slices 
xslice = [-0.5*lx,0.5*lx]; 
yslice = [-0.5*ly,0.5*ly]; 
zslice = [0,Dz/5];

hslice = surf(linspace(0,Dz,100),linspace(-Dx,Dx,100),zeros(100));
rotate(hslice,[-1,0,0],-45)

xd = get(hslice,'XData'); yd = get(hslice,'YData'); zd = get(hslice,'ZData');

figure 
colormap(jet) 
h = slice(z,x,y,abs(Wo),zd,xd,yd); 
h.FaceColor = 'interp'; 
h.EdgeColor = 'none'; 
h.DiffuseStrength = 0.8;

figure(9)
set(gcf, 'color', [0 0 0])
h= slice(z,x,y,abs(Wo),zslice,xslice,yslice);
% set(gca,'Color','k')
zlim([0 Dz])
colormap(mapgreen)
% set(h,'EdgeColor','none',...
% 'FaceColor','interp',...
% 'FaceAlpha','interp')
% alpha('color')

% colormap('hsv')
shading interp
alpha(0.7); %axis off
% set(gcf,'Color',[1 1 1])
% colormap(mapgreen)

daspect([1.5,.1,.1]) 
axis tight 
view(-38.5,16) 
camzoom(1.4) 
camproj perspective
axis off

hold on

%  plotPropagatedRays(rayH11,rayH22);
 hold off
%%
%% distances

distances = [2/2,6/2,inf];

texts     = {'$z_R/2$','$z_R/6$','0'};


FG = W22o;

close(figure(800))
fig800 = figure(800);
ha = tight_subplot(2,3,[.01 .01],[.05 .08],[.12 .09]);


kk = 1;
for jj = distances
  
  index = floor(Nz/jj);
  if index == 0
    index = 1;
  end
  
  gg = W22(:,index,:);
  gg = reshape(gg,[Nx,Nx]);
  axes(ha(4-kk))
  plotOpticalField(x/InitialWaist,x/InitialWaist,abs(gg).^1.5,mapgreen,'');
  title(['$z$ = ', texts{kk}],'Interpreter','latex')

  gg1 = W22o(:,index,:);
  gg1 = reshape(gg1,[Nx,Nx]);
  axes(ha(7 - kk))
  plotOpticalField(x/InitialWaist,x/InitialWaist,abs(gg1).^1.5,mapgreen,'');
  title(['$z$ = ', texts{kk}],'Interpreter','latex')

  kk = kk+1;




  
end

sgtitle('Propagation of $\psi_{n,m}^{(1,1)}$','interpreter','latex')


ha(1).XAxis.Visible = 'off';
ha(2).XAxis.Visible = 'off';
ha(3).XAxis.Visible = 'off';
ha(2).YAxis.Visible = 'off';
ha(5).YAxis.Visible = 'off';
ha(3).YAxisLocation = 'right';
ha(6).YAxisLocation = 'right';
% saveas(gcf,['HH11Propagation.png'])
export_fig('HH11Propagation-gauss','-png','-transparent')


%%

distances = [25/12,10/2,inf];

texts     = {'6/25','1/10','0'};


FG = W22o;

close(figure(800))
fig800 = figure(800);
ha = tight_subplot(2,3,[.01 .01],[.05 .1],[.1 .07]);


kk = 1;
for jj = distances
  
  index = floor(Nz/jj);
  if index == 0
    index = 1;
  end
  

  gg1 = W11o(:,index,:);
  gg1 = reshape(gg1,[Nx,Nx]);
  figure(1000)
  plotOpticalField(x,x,abs(gg1).^1.5,mapgreen,'');
  plotRays(rayH21(index),'y',1)
  title(['$z$ = ', texts{kk}],'Interpreter','latex')

  kk = kk+1;

  saveas(gcf,['HH22Propagation',num2str(jj),'.png'])

  
  gg1 = W12o(:,index,:);
  gg1 = reshape(gg1,[Nx,Nx]);
  figure(1000)
  plotOpticalField(x,x,abs(gg1).^1.5,mapgreen,'');
  plotRays(rayH21(index),'y',1)
  title(['$z$ = ', texts{kk}],'Interpreter','latex')

  kk = kk+1;

  saveas(gcf,['HH12Propagation',num2str(jj),'.png'])

  
  

  
end

%%

distances = [1,3,inf];

texts     = {'1/5','1/10','0'};

kk = 1;
for jj = distances
  
  index = floor(Nz/jj);
  if index == 0
    index = 1;
  end
  
  gg1 = W(:,index,:);
  gg1 = reshape(gg1,[Nx,Nx]);
  figure(1000)
  plotOpticalField(x/InitialWaist,x/InitialWaist,abs(gg1).^1.5,mapgreen,'');
%   plotRaysSquare(rayH12(index),'m',InitialWaist)
%   plotRaysSquare(rayH21(index),'y',InitialWaist)
%   plotRaysSquare(rayH11(index),'r',InitialWaist)
%   plotRaysSquare(rayH22(index),'c',InitialWaist)
  title(['$z$ = ', texts{kk}],'Interpreter','latex')

  kk = kk+1 

%   saveas(gcf,['HermitePropagation',num2str(jj),'.png'])
export_fig(['ImagesThesis/HermitePropagation',num2str(jj)],'-png','-transparent')


  
end




%%


% disz = 1,...floor(Nz/3);
% gg = W(:,disz,:);
% %gg = g(:);
% gg = reshape(gg,[Nx,Nx]);
% figure(8)
% plotOpticalField(x/InitialWaist,x/InitialWaist,abs(gg).^1.5,mapgreen,'');
% title(['$\zeta$ = ', num2str(z(disz)/RayleighDistance)],'Interpreter','latex')
% 
% 
% plotRays(rayH11(disz),'r',InitialWaist)
% plotRays(rayH21(disz),'y',InitialWaist)
% plotRays(rayH12(disz),'m',InitialWaist)                                         
% plotRays(rayH22(disz),'c',InitialWaist)
% 
% 
% saveas(gcf,'HGBOzR0r.png')
% 
% % rayH11s = rayH11;
% % rayH21s = rayH21;
% % rayH12s = rayH12;
% % rayH22s = rayH22;


%%
% 
% figure(100)
% imagesc(abs(gx))
%% last figure of paper


