%--------------propagacion paraxial por espectro angular------------------%
mapgreen = AdvancedColormap('kg',256,[0 255]/255);  %color del haz

%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Numero de puntos
n=-N/2:N/2-1;               % Vector indicial igualmente espaciado

%muestreo de los vectores espaciales
Dx=20;          % Tamano de ventana para el obstáculo centrado en 0
dx=Dx/N;        % Resolución
x=n*dx;         % Vector x (plano de campo a propapar)
[X]=meshgrid(x);
%muestreo del vector de frecuencia
Du=1/dx;        % Ventana de Fourier
du=1/Dx;        % Resolución
u=n*du;         % Vector
[U]=meshgrid(u);
%vectores kx,ky
kx=2*pi*u;
[Kx]=meshgrid(kx);
%-------------------------------------------------------------------------% 

%--------------Hermite Gauss normalizado(sin unidades físicas)-----------&

%--------------------Hermite Gauss sin ostrucción------------------------%
%indices del Laguerre Gauss
nu=10;mu=9;

%Las dos soluciones del Hermite Gauss para x e y
[HGx,NHGx]= hermitegauss(nu,x);
[HGy,NHGy]= hermitegauss(mu,x);


%funciones Hankel 1 
H1x=HGx+i*NHGx;
H1y=HGy+i*NHGy;
%funciones Hankel 2
H2x=HGx-i*NHGx;
H2y=HGy-i*NHGy;

figure(1)
plot(abs(NHGy))
