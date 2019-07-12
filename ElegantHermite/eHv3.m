%eHermite analiticos
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=50;mu=50;
% Parametros fisicos [micras]
wo=2;
lamb=0.6328;    %micras
k=2*pi/lamb;
zl=k*wo^2/2;
qo=i*zl;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl;%/5; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^5;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=0:dz:Dz;       % Vector z de propagacion
% calculando el tamaño de zp
sizezp=size(z);
% Cintura de gaussiana a la distancia del tamaño de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2);
% Cintura de Laguerre a la distancia del tamaño de la ventana en z
niu=max(nu,mu);
sigmaH=ws*sqrt((2*niu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tamaño de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(sigmaH)*2;%*1.6/(1.97)^2;  % Tamaño de la ventana
dx=Dx/N;            % Resolución
x=n*dx;             % Vector x (plano de campo a propapar)
y=x;
[X]=meshgrid(x);
x1=x(N/2:N);    % vector de 0 a Dx
% Muestreo del vector de frecuencia
Du=1/dx;        % Tamaño de la ventana de Fourier
du=1/Dx;        % Resolución
u=n*du;         % Vector
[U]=meshgrid(u);
% Vectores kx,ky
kx=2*pi*u;
[Kx]=meshgrid(kx); 
% Para dibujar circulos 
ang=0:0.01:2*pi;
%----------------Cintura de elegant Hermite en z=0------------------------%
%--------curvas con z=constante=0
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=x/wo;
% yi=x/wo;
% zi=0;
% Parametros de eHermite
% qz=zi+i*zl;    C=1;
% An=C*(qo./qz).^((nu+1/2));
% Am=C*(qo./qz).^((mu+1/2));
% An=1;Am=1;
% alphax=sqrt(i*(k)./(2*qz)).*xi;
% alphay=sqrt(i*(k)./(2*qz)).*yi;
alphax=xi;
alphay=xi;
[eHGzx,eNHGzx]=ehermite(nu,(alphax));
[eHGzy,eNHGzy]=ehermite(mu,(alphay));
% eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
% eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy; 
%Funciones Hankel de elegant Hermite%
%Funcion Hankel con z=cte
H1z=(eHGzy')*(eHGzx)+i*(eNHGzy')*(eNHGzx);
H2z=(eHGzy')*(eHGzx)-i*(eNHGzy')*(eNHGzx);
eH=(eHGzy')*(eHGzx);
figure(1)
imagesc(abs(eH))
figure(2)
imagesc(angle(eH))