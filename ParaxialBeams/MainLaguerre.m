

%adding guassian parameters folders in subfolder
here = mfilename('fullpath');
[path, ~, ~] = fileparts(here);
addpath(genpath(path));

%Color of beam
mapgreen = AdvancedColormap('kgg',256,[0 70 255]/255);  %color of beam

%-------------------- indices of Laguerre Gaussian Beams -----------------%
nu     = 8;
mu     = 2;

% Physical parameters [microns]
wo     = 1000;
lamb   = 0.6328;
k      = 2*pi/lamb;
zo     = k*wo^2/2;

%------------------------ sampling of vectors ----------------------------%
%First we estimate samplig in z-direction with propagation distance 
% z-direction
Dz     = 2*zo;                    % z-window (propagation distance)
Nz     = 2^6;                     % number of points in z-direction
dz     = Dz/Nz;                   % Resolution in z
z      = 0:dz:Dz;                 % z-vector z of propagation 

% waist of Gaussian until z-propagation
ws     = waistPhysicalGaussianBeam(Dz,wo,zo);
% waist of Gaussian Laguerre Beam until z-propagation
sigmaL = ws*sqrt(2*(2*nu+mu+1));

%Second we estimage sampling in x,y-direction in terms of waist of Guassian
%Laguerre Beam

% y,x-direction
N      =  2^9;                    % Number of points in x,y axis
n      = -N/2+.05:N/2-1+.05;      % vector with N-points with resolution 1
Dx     = (sigmaL);                % Size of window 
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

%Polar coordidates for Laguerre Gauss Beam
[TH,Rho] = cart2pol(X,X');

%% ----------------------- Laguerre Gauss in z = 0 --------------------- %%
zinit  = 0;
% Solutions of Laguerre
LGB    =  laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rho,TH,zinit);

XGB    = xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rho,TH,zinit);

% Hankel Functions
H1     = LGB+1i*XGB;
H2     = LGB-1i*XGB;

% Function to propagate 
g      = LGB;
% Plot of Function
figure(1)
pcolor(x/(sqrt(2)*wo), x/(sqrt(2)*wo), abs(g))
axis square
shading flat
colormap(mapgreen)
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
% Max Peak
pxy     = max(max(g));

%------------------ Obstruction on Lagurre in z = 0 ----------------------%
sigmaLo = wo*sqrt(2*(2*nu+mu+1));   % waist of Laguerrre in z = 0
lo      = sigmaLo/5;                % size of obstruction in terms of waist of Laguerre
xt      = 0;                        % traslation of obstruction in x-axis
yt      = 0;                        % traslation of onstruction in y-axis
[~,rho]  = cart2pol(X-xt,X'-yt);    % Convert this in polar coordinates
obo     = double(rho<=lo);          % Create Obstruction   
clear ro                            % Clean Matrix of polar coordinates
% Applying obstructi obstruction
g       = g.*(1-obo);
%Ploting Laguerre
figure(2)
pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)


%% ----------------------- Ray tracing (rx,z=0)  ----------------------- %%

pn  =  30;                          % Number of rays
ray = ([]);                         % Structure of rays

for jj=1:pn
    %vectors of zeros for each x,y-component ray for H1 and H2
    ray(jj).rxH1 = zeros(1,length(z));
    ray(jj).rxH2 = zeros(1,length(z));
    ray(jj).ryH1 = zeros(1,length(z));
    ray(jj).ryH2 = zeros(1,length(z));
end


for jj = 1:pn
    % Cartersian coordinates of point in circunference of obstruction
    xi              = xt+lo*cos(jj*(2*pi)/(pn)); 
    yi              = yt+lo*sin(jj*(2*pi)/(pn));
    % saving cartesian coordinates of point in circunference in each ray
    ray(jj).xc      = xi; 
    ray(jj).yc      = yi;
    ray(jj).rxH1(1) = xi; 
    ray(jj).ryH1(1) = yi;
    ray(jj).rxH2(1) = xi; 
    ray(jj).ryH2(1) = yi;
    
    % H1 with y = cte, z = cte
    Rhop = sqrt(x.^2+yi.^2);
    Thp  = atan2(yi,x);
    H1x  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1))...
         + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1));
    fx   = unwrap(angle(H1x));
    
    % H1 with x = cte, z = cte
    Rhop = sqrt(xi.^2+y.^2);
    Thp  = atan2(y,xi);
    H1y  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1))...
         + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1));
    fy   = unwrap(angle(H1y));
    
    % H1 with x = cte, y = cte
    Rhop = sqrt(xi.^2+yi.^2);
    Thp  = atan2(yi,xi);
    H1z  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z)...
         + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z);
    fz   = unwrap(angle(H1z));  
            
    % Calculating slopes of H1
    [ray(jj).mzxH1,ray(jj).mzyH1,ray(jj).mxyH1] = gradientxyz(fx,fy,fz,k,dx,dx,dz,xi,yi,z(1));

    % H2 with y = cte, z = cte
    Rhop = sqrt(x.^2+yi.^2);
    Thp  = atan2(yi,x);
    H2x  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1))...
         - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1));
    fx   = unwrap(angle(H2x));
    
    % H2 with x = cte, z = cte
    Rhop = sqrt(xi.^2+y.^2);
    Thp  = atan2(y,xi);
    H2y  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1))...
         - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(1));
    fy   = unwrap(angle(H2y));
    
    % H2 with x = cte, y = cte
    Rhop = sqrt(xi.^2+yi.^2);
    Thp  = atan2(yi,xi);
    H2z  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z)...
         - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z);
    fz   = unwrap(angle(H2z));  
            
    % Calculating slopes of H2
    [ray(jj).mzxH2,ray(jj).mzyH2,ray(jj).mxyH2] = gradientxyz(fx,fy,fz,k,dx,dx,dz,xi,yi,z(1));               
 
end

% Initial Field with rays in this init conditions
figure(3)
pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^2)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)
hold on
for jj=1:pn
    plot(ray(jj).xc/(sqrt(2)*wo),ray(jj).yc/(sqrt(2)*wo)...
        ,'.','MarkerSize',10,'LineWidth',2,'color','r')
end
hold off


%%  ----------------------- Physical Propagation ------------------------ %
% paraxial propagator 
prop = exp(-1i*lamb*dz*(Kx.^2+(Kx').^2)/(4*pi));
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
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^1)
    axis off
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    drawnow 
   hold on
   % propagated points of H1 and H2
    for jj=1:pn
        plot(ray(jj).rxH1(ii-1)/(sqrt(2)*wo),ray(jj).ryH1(ii-1)/(sqrt(2)*wo),'.','MarkerSize',20,'LineWidth',2,'color','r')
        plot(ray(jj).rxH2(ii-1)/(sqrt(2)*wo),ray(jj).ryH2(ii-1)/(sqrt(2)*wo),'.','MarkerSize',20,'LineWidth',2,'color','y')
    end
   hold off
   pause(1)

    %------------------------ Calculating Rays ---------------------------%   
    % Given z we find z+dz and new position x+dx and y+dy with help of
    % slopes

    for jj=1:pn
        % step of ray in z-direction, equation or ray r(z) = mrz*dz+r(z-1)
        % dependece of last point, new vector
        ray(jj).rxH1(ii)   = ray(jj).rxH1(ii-1) + (1/ray(jj).mzxH1)*(z(ii)-z(ii-1));
        ray(jj).ryH1(ii)   = ray(jj).ryH1(ii-1) + (1/ray(jj).mzyH1)*(z(ii)-z(ii-1));
        
        ray(jj).rxH2(ii)   = ray(jj).rxH2(ii-1) + (1/ray(jj).mzxH2)*(z(ii)-z(ii-1));
        ray(jj).ryH2(ii)   = ray(jj).ryH2(ii-1) + (1/ray(jj).mzyH2)*(z(ii)-z(ii-1));
  
        % ------------------ Calculating Slopes ------------------------- %

        % H1 with y = cte, z = cte
        Rhop = sqrt(x.^2+ray(jj).ryH1(ii).^2);
        Thp  = atan2(ray(jj).ryH1(ii),x);
        H1x  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii))...
             + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii));
        fx   = unwrap(angle(H1x));

        % H1 with x = cte, z = cte
        Rhop = sqrt(ray(jj).rxH1(ii) .^2+y.^2);
        Thp  = atan2(y,ray(jj).rxH1(ii));
        H1y  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii))...
             + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii));
        fy   = unwrap(angle(H1y));

        % H1 with x = cte, y = cte
        Rhop = sqrt(ray(jj).rxH1(ii).^2+ray(jj).ryH1(ii).^2);
        Thp  = atan2(ray(jj).ryH1(ii),ray(jj).rxH1(ii));
        H1z  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z)...
             + 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z);
        fz   = unwrap(angle(H1z));  

        % Calculating slopes of H1
        [ray(jj).mzxH1,ray(jj).mzyH1,ray(jj).mxyH1] = gradientxyz(fx,fy,fz,k,dx,dx,dz,ray(jj).rxH1(ii),ray(jj).ryH1(ii),z(ii));

        % H2 with y = cte, z = cte
        Rhop = sqrt(x.^2+ray(jj).ryH2(ii).^2);
        Thp  = atan2(ray(jj).ryH2(ii),x);
        H2x  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii))...
             - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii));
        fx   = unwrap(angle(H2x));

        % H2 with x = cte, z = cte
        Rhop = sqrt(ray(jj).rxH2(ii).^2+y.^2);
        Thp  = atan2(y,ray(jj).rxH2(ii));
        H2y  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii))...
             - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z(ii));
        fy   = unwrap(angle(H2y));

        % H2 with x = cte, y = cte
        Rhop = sqrt(ray(jj).rxH2(ii).^2+ray(jj).ryH2(ii).^2);
        Thp  = atan2(ray(jj).ryH2(ii),ray(jj).rxH2(ii));
        H2z  =     laguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z)...
             - 1i*xLaguerrePhysicalGaussBeam(nu,mu,wo,zo,Rhop,Thp,z);
        fz   = unwrap(angle(H2z));  

        % Calculating slopes of H2
        [ray(jj).mzxH2,ray(jj).mzyH2,ray(jj).mxyH2] = gradientxyz(fx,fy,fz,k,dx,dx,dz,ray(jj).rxH2(ii),ray(jj).ryH2(ii),z(ii));               
       
    end
    %------------------------ End calculating rays -----------------------%   
    %propagating field
    G=fftshift(fft2(g));
    %obtain new propagated field
    g=ifft2(fftshift(G.*prop));
    %saving transversal fields
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);
%     
end

