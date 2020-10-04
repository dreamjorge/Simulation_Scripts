
%%          Script of Laguerre Beam (properties and propagation)
% adding path for classes and functions
addpath ParaxialBeams
addpath ParaxialBeams\Addons
% Selecting green color for beam
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);
% use gpu? 
runGPU = 'no';
%%          Initial parameters of Laguerre Gaussian Beams

l = 11;
p = 0;

% Physical parameters [microns]
InitialWaist        = 10;%179*2;
Wavelength          = 0.6328;
% normalized parameters
% InitialWaist = 1;
% Wavelength   = pi;
PropagationDistance = 0;

% Calculating Laguerre parameters in z=0
LPinZ0 = LaguerreParameters(PropagationDistance...
                           ,InitialWaist...
                           ,Wavelength...
                           ,l...
                           ,p);
                     
%%                        Sampling of vectors 
% Estimate sampling in z-direction with propagation distance 
% z-direction
Dz = LPinZ0.RayleighDistance; % z-window (propagation distance)
Nz = 2^9;                     % number of points in z-direction
dz = Dz/Nz;                   % Resolution in z
z  = 0:dz:Dz;                 % z-vector z of propagation 

% Calculating Laguerre parameters in z distance vector
% LPinZ = LaguerreParameters(z...
%                           ,InitialWaist...
%                           ,Wavelength...
%                           ,l...
%                           ,p);
                        
% Calculating Laguerre parameters in z distance vector copying object in z = 0
LPinZ  = copy(LPinZ0);
LPinZ.zCoordinate = z;                         
                        
fig1          = figure(1);
fig1.Position = [314 300 1097 479];
plotLaguerreParameters(LPinZ);

% Estimate sampling in x,y-direction in terms of waist of max  
% Laguerre Gauss Beam waist until max z-propagation

% MaxLaguerreWaist = LaguerreParameters.waistFunction(z(end)...
%                                                    ,InitialWaist...
%                                                    ,LPinZ0.RayleighDistance...
%                                                    ,l...
%                                                    ,p);
 % waist in max z, from LaguereParametersz0                                               
MaxLaguerreWaist = LPinZ.LaguerreWaist(end);

N   = 2^10;                  % Number of points in x,y axis
n   = -N/2+.05:N/2-1+.05;   % vector with N-points with resolution 1
Dx  = 2.2*MaxLaguerreWaist; % Size of window 
dx  = Dx/N;                 % Resolution
x   = n*dx;                 % Vector with dimentions
y   = x;                    % same for y
[X] = meshgrid(x,y);        % Matrix for x,y


% Transformation of coordinates
[dth,dr]                      = cart2pol(dx,dx');
[thetaCoordinate,rCoordinate] = cart2pol(x,x);
[ThetaCoordinate,RCoordinate] = cart2pol(X,X');

r  = 0:N-1;
Dr = sqrt(2)*Dx;
dr = Dr/N;
r  = r*dr;

th   = -N/2:N/2-1;
Dth  = atan2(Dx,Dx);
dth  = Dth/N;
th   = th*dth;

% Diferential vector in cartesian coordinates
%difr = [dx,dx,dz];
difr = [dr,dth,dz];
% Estimate vectors of frequency for Fourier Transforms associated
% with x,y
Du   = 1/dx;         % Size of window 
du   = 1/Dx;         % Resolution
u    = n*du;         % freq vector with dimentions
[U]  = meshgrid(u);  % freq Matrix
%kx,ky vectors
kx   = 2*pi*u;       % angular freq vector
[Kx] = meshgrid(kx); % angular freq matrix


if strcmp(runGPU,'yes')
%Parallel Computing Toolbox
  x               = gpuArray(x);
  X               = gpuArray(X);
  u               = gpuArray(u);
  U               = gpuArray(U);
  kx              = gpuArray(kx);
  rCoordinate     = gpuArray(rCoordinate);
  thetaCoordinate = gpuArray(thetaCoordinate);
  RCoordinate     = gpuArray(RCoordinate);
  ThetaCoordinate = gpuArray(ThetaCoordinate);
end
%% ----------------------- Laguerre Gauss in z = 0 --------------------- %%

% With laguerre parameters calculated, it estimates Laguerre Gauss Beam
LG = LaguerreBeam(RCoordinate,ThetaCoordinate,LPinZ0);

% Plot of Field
figure(2)
subplot(1,2,1)
plotOpticalField(x,x,abs(LG.OpticalFieldLaguerre),mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);

% Optic Field to propagate 
g   = LG.OpticalFieldLaguerre;
% Max Peak
pxy = max(max(g));

%% ----------------- Obstruction on Lagurre in z = 0 ------------------- %%
% Initial Waist of Laguerre Beam
% 
% lo      = (LPinZ0.LaguerreWaist)/4.3;  % size of obstruction in terms of waist of Laguerre
% xt      = LPinZ0.LaguerreWaist/5.8;                           % traslation of obstruction in x-axis
% yt      = LPinZ0.LaguerreWaist/5.8;                           % traslation of onstruction in y-axis
lo      = (LPinZ0.LaguerreWaist)*.195;  % size of obstruction in terms of waist of Laguerre
xt     = 0;                           % traslation of obstruction in x-axis
yt     = 0;         

[~,rho] = cart2pol(X-xt,X'-yt);        % Convert this in polar coordinates
obo     = double(rho<=lo);             % Create Obstruction   
clear rho      
% Clean Matrix of polar coordinates
% Applying obstruction in optic field
g       = g.*(1-obo);
%Ploting Laguerre with obstruction
figure(3)
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);
plotCircle(xt,yt,lo);


%% ----------------------- Ray tracing (rx,z=0)  ----------------------- %%

TotalRays = 4;          % Number of rays

rayH1(Nz) = CylindricalRay();
rayH2(Nz) = CylindricalRay();

for point_index = 1 : TotalRays
    % Cartersian coordinates of point in circunference of obstruction
    xi = xt + lo*cos(point_index*(2*pi)/(TotalRays))+.001; 
    yi = yt + lo*sin(point_index*(2*pi)/(TotalRays))+.001;
    zi = 0;
    % assign coordinate to Optical Rays in z = 0, i.e index_z = 1  
    [rayH1(1)] = assignCoordinates2CylindricalRay(xi,yi,zi,rayH1(1),point_index,1);
    [rayH2(1)] = assignCoordinates2CylindricalRay(xi,yi,zi,rayH2(1),point_index,2);
    
end

% Initial Field with rays in this init conditions
figure(3) 
plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
plotRays(rayH1(1),'r',1)


%%                         Physical Propagation
% paraxial propagator 
prop    = paraxialPropagator(Kx,Kx',LPinZ0.k,dz);

figure(5)
imagesc(u,u,(angle(prop)))
title('Propagator')

% Matrix for save transversal fields
gx      = zeros(N,length(z)); 
gy      = zeros(N,length(z));

if strcmp(runGPU,'yes')
  gx      = gpuArray(gx); 
  gy      = gpuArray(gy); 
end

% Save field in z = 0 
gx(:,1) = g(N/2+1,:);
gy(:,1) = g(:,N/2+1);


for z_index = 1:length(z)-1 % propagation with respect to z
%% loop of each component in z

  % propagation of Optical Field 
  g = propagateOpticalField(g,prop);
  %saving transversal fields
  gx(:,z_index+1) = g(N/2+1,:);
  gy(:,z_index+1) = g(:,N/2+1);  
    
  %%                    Calculating Rays

  % propagation distance 
  zi     = z(z_index);
  % calculating Laguerre Parameters in zi
  LPinZi = LaguerreParameters(zi,InitialWaist,Wavelength,l,p);   

  % propagate all rays of H1
  [rayH1(z_index+1)] = getPropagateCylindricalRays(rayH1(z_index),...
                                                   TotalRays,...
                                                   r,th,...
                                                   difr,...
                                                   LPinZi,...
                                                   LPinZ...
                                                   ); 
 

  % propagate all rays of H2
  [rayH2(z_index+1)] = getPropagateCylindricalRays(rayH2(z_index),...
                                                   TotalRays,...
                                                   r,th,...
                                                   difr,...
                                                   LPinZi,...
                                                   LPinZ...
                                                   );
                         

  %                             End calculating rays 
  %%

  fig = figure(6);
  fig.Position = [-1349 147 813 733];
  plotOpticalField(x,x,abs(g).^2,mapgreen,'microns');
  title(['z = ', num2str(z_index), ' of ', num2str(Nz)])
  drawnow 
  hold on

  % Plot propagated points of H1 and H2 in iteration before
  plotRays(rayH1(z_index+1),'r')
  plotRays(rayH2(z_index+1),'y')

  pause(0.1)

end
%%
figure(7)
plotRaysPropagated(rayH1,rayH2,Nz);
%%
figure(8)
imagesc(z,x,abs(gx))
hold on
for z_index = 1:length(z)-1
  scatter(rayH2(z_index).zCoordinate,rayH2(z_index).xCoordinate,10,'filled','MarkerFaceColor',[1 0 0])
  scatter(rayH1(z_index).zCoordinate,rayH1(z_index).xCoordinate,10,'filled','MarkerFaceColor',[0 1 0])
end
hold off

figure(9)

for z_index = 1:length(z)-1
  
  radial1(z_index) = rayH2(z_index).rCoordinate(1);
  radial2(z_index) = rayH2(z_index).rCoordinate(2);
end
hold on
plot(radial1)
plot(radial2)
hold off