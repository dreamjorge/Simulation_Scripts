%% physical units in microns

InitialWaist  = 0.00005;
Wavelength    = 6.33e-5;
K             = 2*pi/Wavelength;
%% elegant normalization 
InitialWaist  = 1;
Wavelength    = pi;
K             = 2*pi/Wavelength;
%% z - direction vector

N      =  2^8;                    % Number of points in z direction
n      = (0:N-1);                 % N
Dz     = 2;                       % Size of window 
dz     = Dz/N;                    % Resolution
z      = n*dz;                    % Vector

%% Estimate Gaussian Parameters with this z-vector
GP = GaussianParameters(z,InitialWaist,Wavelength);
%  plot waist
figure(1)
plot(GP.PropagationDistance,GP.Waist) 
hold on
plot(GP.PropagationDistance,-GP.Waist)
hold off

%% vector in x,y using waist of beam in z max
N      =  2^9;                    % Number of points in x,y axis
n      = (-N/2:N/2-1);            % vector with N-points with resolution 1
Dx     = 4*GP.Waist(end);         % Size of window 
dx     = Dx/N;                    % Resolution
x      = n*dx;                    % Vector
y      = x;
[X,Y]  = meshgrid(x,y);
[TH,R] = cart2pol(X,Y);
% u,v vector for spectrum Fourier
Du     = 1/dx;                    % Size of window 
du     = 1/Dx;                    % Resolution
u      = n*du;                    % Vector
[U]    = meshgrid(u);
%kx,ky vectors
kx     = 2*pi*u;
[Kx]   = meshgrid(kx);


%% Gaussian in z = 0 and propagation
GB              = GaussianBeam(X,Y,0,InitialWaist,Wavelength);
amplitudefactor = max(GBA.OpticalField(:));
% this field we'll propagate
g               = GB.OpticalField;

% plot this field
figure(2)
imagesc(x,y,abs(g).^2)
caxis([0 amplitudfactor])
hold on
%plot circle of waist
r     = GP.Waist(1);
th    = 0:pi/50:2*pi;
xt    = 0;
yt    = 0;
xunit = r * cos(th) + xt;
yunit = r * sin(th) + yt;
plot(xunit, yunit,  'LineWidth',2,'Color','r');
hold off
axis square

% propagator
prop=exp(1i*Wavelength*dz*(Kx.^2+(Kx').^2)/(4*pi));
figure(3)
imagesc(u,u,(angle(prop)))
title('Propagator')

%%
for ii = 2:length(z) % propagation with respect to z
  
    %% field before propagation
    figure(4)
    subplot(2,2,1)
    imagesc(x,y,abs(g).^2)
    title('Amplitude Angular Spectrum Method')
    caxis([0 amplitudfactor])
    hold on
    %plot circle of waist
    r     = GP.Waist(ii);
    th    = 0:pi/50:2*pi;
    xt    = 0;
    yt    = 0;
    xunit = r * cos(th) + xt;
    yunit = r * sin(th) + yt;
    plot(xunit, yunit,  'LineWidth',2,'Color','r');
    hold off
    axis square
    
    subplot(2,2,2)
    gp = GaussianBeam(X,Y,z(ii),InitialWaist,Wavelength);
    imagesc(x,y,abs(gp.OpticalField).^2)
    title( 'Amplitude Analytic Function')
    caxis([0 amplitudfactor])
    hold on
    plot(xunit, yunit,  'LineWidth',2,'Color','r');
    hold off
    axis square
    
    subplot(2,2,3)
    imagesc((angle(g)))
    title('Phase of Propagated Field')
    axis square
    
    subplot(2,2,4)
    imagesc((angle(gp.OpticalField)))
    title('Phase of Analytic Field')
    
    axis square
    drawnow
    
    % propagation
    % it's needed correction in phase of FFT
    G =fftshift(fft2(ifftshift(g)));...*exp(-1j.*(x(1)).*(Kx));...*exp(-1j.*(x(1)).*(Kx'));
    %obtain new propagated field
    g=fftshift(ifft2(ifftshift(G.*prop)));...*exp(-1j.*(u(1)).*(X));...*exp(-1j.*(u(1)).*(X'));
  
    pause(.5)
end

