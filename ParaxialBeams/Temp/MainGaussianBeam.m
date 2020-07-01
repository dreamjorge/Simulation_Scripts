wo      = 0.1191;             % Initial normalized standard waist of Gaussian Beam.
sigmao  = 2*wo;          % Initial normalized deviation of Gaussian Beam.
weo     = wo/sqrt(2);    % Initial normalized elegant waist of Gaussian Beam.
% Normalized distance of Rayleigh
so      = 40;


% Quantities for generate vector in s direction
Ns      = 2^9;            % Number of points of vector
ns      = -Ns/2:Ns/2-1;   % Index vector with 1 of resolution
timesso = 12;             % Number of times so 
Ds      = timesso*so;     % Size of vector's window
ds      = Ds/Ns;          % Resolution of vector
s       = ns.*ds;         % Vector


% y,x-direction
N      =  2^9;                    % Number of points in x,y axis
n      = -N/2+.05:N/2-1+.05;      % vector with N-points with resolution 1
Dx     = 4*wo;                    % Size of window 
dx     = Dx/N;                    % Resolution
x      = n*dx;                    % Vector
y      = x;
[X,Y]  = meshgrid(x,y);

GB = GaussianBeam(X,Y,0,wo,0.065);
figure(1)
imagesc(x,y,abs(GB.OpticalField))
