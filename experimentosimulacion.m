load('C:\Users\jorge\Dropbox\AJwork\Imagenes lateral Laguerre Experimental Procesado\Laguerre Obstruccion fuera de eje\1.mat')
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
g  = LG.OpticalFieldLaguerre;
% Plot of Field
figure(2)
subplot(1,2,1)
plotOpticalField(x,x,abs(g),mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);

% Max Peak
pxy = max(max(g));
obradius = LG.InitialWaist*.23; 

figure(2)
subplot(1,2,2)
imagesc(abs(im))
axis square
waistpixels = 128;
xc = 844/2-1;
yc = 844/2;
plotCircle(xc,yc,waistpixels);
xt = xc + waistpixels/8;
yt = yc + 0; 
obstructionradius = waistpixels*.23;
plotCircle(xt,yt,obstructionradius);



lo      = (LPinZ0.LaguerreWaist)*.2;  % size of obstruction in terms of waist of Laguerre
xto     = (LPinZ0.LaguerreWaist)/8;                           % traslation of obstruction in x-axis
yto     = 0;                           % traslation of onstruction in y-axis
[~,rho] = cart2pol(X-xto,X'-yto);        % Convert this in polar coordinates
obo     = double(rho<=lo);             % Create Obstruction   
clear rho      
% Clean Matrix of polar coordinates
% Applying obstruction in optic field
g1       = g.*(1-obo);
figure(2)
subplot(1,2,1)
plotOpticalField(x,x,abs(g1),mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);
plotCircle(xto,yto,lo);

figure(3)
subplot(1,2,2)
imagesc(abs(cdata));
axis square
waistpixels = 129;
xc = 1008/2 - 25;
yc = 1008/2 + 39;
plotCircle(xc,yc,waistpixels);
xt = xc + waistpixels*.76;
yt = yc + 0; 
obstructionradius = waistpixels*.23;
plotCircle(xt,yt,obstructionradius);
lo      = (LPinZ0.LaguerreWaist)*.22;  % size of obstruction in terms of waist of Laguerre
xto     = (LPinZ0.LaguerreWaist)*.78;                           % traslation of obstruction in x-axis
yto     = 0;                           % traslation of onstruction in y-axis
[~,rho] = cart2pol(X-xto,X'-yto);        % Convert this in polar coordinates
obo     = double(rho<=lo);             % Create Obstruction   
clear rho      
% Clean Matrix of polar coordinates
% Applying obstruction in optic field
g2       = g.*(1-obo);
subplot(1,2,1)
plotOpticalField(x,x,abs(g2),mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);
plotCircle(xto,yto,lo);

% central
lo      = (LPinZ0.LaguerreWaist)*.195;  % size of obstruction in terms of waist of Laguerre
xto     = 0;                           % traslation of obstruction in x-axis
yto     = 0;                           % traslation of onstruction in y-axis
[~,rho] = cart2pol(X-xto,X'-yto);        % Convert this in polar coordinates
obo     = double(rho<=lo);             % Create Obstruction   
clear rho      
% Clean Matrix of polar coordinates
% Applying obstruction in optic field
g2       = g.*(1-obo);
subplot(1,2,1)
plotOpticalField(x,x,abs(g2).^2,mapgreen,'microns');
plotCircle(0,0,LPinZ0.LaguerreWaist);
plotCircle(xto,yto,lo);