%% add path for classes and functions

addpath ParaxialBeams
addpath ParaxialBeams\Addons
mapgreen = AdvancedColormap('kgg',256,[0 100 255]/255);  %color of beam

% Physical parameters [microns]
InitialWaist           = 100;%179*2;
Wavelength             = 0.6328;

% Normalized units
% InitialWaist         = 1;
% Wavelength           = pi;

% Calculating parameters of Gaussian Beam in z=0
GaussianBeamParameters = GaussianParameters(0,InitialWaist,Wavelength);

% %------------------------ sampling of vectors ----------------------------%
%First we estimate sampling in z-direction with propagation distance 
% z-direction
Dz = 2*GaussianBeamParameters.RayleighDistance;    % z-window (propagation distance)
Nz = 2^8;                                        % number of points in z-direction
dz = Dz/(Nz-1);                                  % Resolution in z
nz = 0:(Nz-1);                                   % Index vector
z  = nz*dz;                                      % z-vector z of propagation 

% Calculating
GaussianBeamParameters = GaussianParameters(z,InitialWaist,Wavelength);

figure(1)
p1=plot(z,GaussianBeamParameters.Waist,'Color','red');
hold on
p2=plot(z,-GaussianBeamParameters.Waist,'Color','red');
p3=plot(z,z*tan(GaussianBeamParameters.DivergenceAngle),'Color','blue');
p4=plot(z,-z*tan(GaussianBeamParameters.DivergenceAngle),'Color','blue');
xlabel('Distance of Propagation [microns]')
ylabel('[microns]')
title('Parameters of Gaussian Beam')
legend([p1,p3],{'Waist of Gaussian Beam',...
               ['Angle of Divergence = ',num2str(rad2deg(GaussianBeamParameters.DivergenceAngle)),'°']})
hold off
%After sampling z vector, we estimage sampling in x,y-direction in terms of
%waist of max waist Gauss Beam until max z-propagation

MaxWaist = GaussianParameters.waistFunction(z(end),InitialWaist,GaussianBeamParameters.RayleighDistance);

N     = 2^10;                  % Number of points in x,y axis
n     = -N/2+.05:N/2-1+.05;    % vector with N-points with resolution 1
Dx    = (1.2)*2*MaxWaist;     % Size of window 
dx    = Dx/N;                  % Resolution
x     = n*dx;                  % Vector with dimentions
y     = x;
[X]   = meshgrid(x,y);
[~,R] = cart2pol(X,X');

%Last we estimate vectors of frequency for Fourier Transforms associated
%with x,y

Du   = 1/dx;                 % Size of window 
du   = 1/Dx;                 % Resolution
u    = n*du;                 % Vector with dimentions
[U]  = meshgrid(u);
%kx,ky vectors
kx   = 2*pi*u;
[Kx] = meshgrid(kx);
%% ------------------------ Gaussian Beam in z = 0 ---------------------- %%


GaussianBeamParameters = GaussianParameters(0,InitialWaist,Wavelength);
GB  = GaussianBeam(R,GaussianBeamParameters);

% Plot of Field
figure(2)
pcolor(x, x, abs(GB.OpticalField).^2)
axis square
shading flat
colormap(mapgreen)
plotCircle(0,0,GaussianBeamParameters.InitialWaist);
xlabel('$x \left[ microns \right]$','Interpreter','latex') 
ylabel('$y \left[ microns \right]$','Interpreter','latex') 

% Optic Field to propagate 
g   = GB.OpticalField;


% Matrix for save transversal fields
gx      = zeros(N,length(z)); 
gy      = zeros(N,length(z));
% Save field in z = 0 
gx(:,1) = g(N/2+1,:);
gy(:,1) = g(:,N/2+1);


for ii = 2:length(z) % propagation with respect to z
    
    fig = figure(6);
    imagesc(x,x,abs(g).^2)
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    drawnow 
    title(['z = ',num2str(z(ii)/GaussianBeamParameters.RayleighDistance)])
    plotCircle(0,0,GB.Waist);
    pause(0.5)

    %propagated theoric field
    GB = GaussianBeam(X,X',z(ii),InitialWaist,Wavelength);...*exp(-1j.*(u(1)).*(X))*exp(-1j.*(u(1)).*(X'));
    GB.Waist
    g = GB.OpticalField;
    %saving transversal fields
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);
%     
end