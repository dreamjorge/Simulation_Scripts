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
nu=17;mu=18;

%Las dos soluciones del Hermite Gauss para x e y
[HGx,NHGx]= hermitegauss(nu,x);
[HGy,NHGy]= hermitegauss(mu,x);

%funciones Hankel 1 en x
H1x=HGx+i*NHGx;
H1y=HGx+i*NHGx;
%funciones Hankel 2 en y
H2x=HGx-i*NHGx;
H2y=HGx-i*NHGx;

%Hermite xy
HG=(HGy')*HGx;
NHG=(NHGy')*NHGy;

% funcion hankel 1 xy
H1xy=(H1y')*H1x;
% funcion hankel 2 xy
H2xy=(H2y')*H2x;

%graficando el Hermite
figure(1)
imagesc(x,x,abs(HG))
title(['Hermite Gauss con n=',num2str(nu),' m=',num2str(mu)])
axis square
colormap(mapgreen)

% Graficando cortes transversales x
figure(2)
subplot(2,1,1);
plot(x,abs(HG(:,N/2)).^2)
title(['Hermite transversal x con n=',num2str(nu),' m=',num2str(mu)])
subplot(2,1,2);
plot(x,abs(NHG(:,N/2)).^2,'r')
title(['NHermite transversal x con n=',num2str(nu),' m=',num2str(mu)])

%-------------------------------------------------------------------------%
%-------------Pendientes a los frente de ondas en la cintura--------------%
%----------------------del Elegant Hermite en z=0-------------------------%
%--------------------funciones analiticas de eH para z=0------------------%
%-------------------------------------------------------------------------%
%--------curvas con z=constante=0
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=x;
yi=x;
zi=0;
% Parametros de eHermite
wo=1;
wzi=wo.*sqrt(1+(zi./zl).^2);
alphax=sqrt(2)*x./wzi;
[eHGzx,eNHGzx]=hermitegauss(nu,(alphax));
[eHGzy,eNHGzy]=hermitegauss(mu,(alphax));
%Funciones Hankel de elegant Hermite%
% Funciones Hankel 1 en x
H1zx=eHGzx+i*eNHGzx;
H1zy=eHGzy+i*eNHGzy;
% Funciones Hankel 2 en y
H2zx=eHGzx-i*eNHGzx;
H2zy=eHGzy-i*eNHGzy;
%-------curva para x,y=constante y todo z
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=sigmaeHx(1);
yi=sigmaeHy(1);
zi=z;
% Parametros de eHermite
qz=zi-i*zl;    C=1;
An=C*((k/zl)^(1/4))*(-i*zl./qz).^((nu+1)/2);
Am=C*((k/zl)^(1/4))*(-i*zl./qz).^((mu+1)/2);
alphax=sqrt(-i*(k*xi.^2)./(2*qz));
alphay=sqrt(-i*(k*xi.^2)./(2*qz));
[eHGxz,eNHGxz]=hermitegauss(nu,(alphax));
[eHGyz,eNHGyz]=hermitegauss(mu,(alphay));
% Funciones Hankel de elegant Hermite%
% Funciones Hankel 1 en x
H1xz=eHGxz+i*eNHGxz;
H1yz=eHGyz+i*eNHGyz;
% Funciones Hankel 2 en y
H2xz=eHGxz-i*eNHGxz;
H2yz=eHGxz-i*eNHGxz;
% Calculando la pendiente con las derivadas de las fases
% fase z constante y todo x,y
phzx=unwrap(angle(H2zx));
phzy=unwrap(angle(H2zy));
% fase x,y constante y todo z
phxz=unwrap(angle(H2xz));
phyz=unwrap(angle(H2yz));
%derivadas para x,y constantes
wavefrontxz=gradient(phxz)/dz+k;
wavefrontyz=gradient(phyz)/dz+k;
%derivadas para z constante
wavefrontzx=gradient(phzx)/dx;
wavefrontzy=gradient(phzy)/dx;
%calculo de la pendiete para x,y en z=0
mxo=wavefrontxz(floor(z(1)/dz+1))/wavefrontzx(N/2+1+floor(xi/dx));
myo=wavefrontyz(floor(z(1)/dz+1))/wavefrontzy(N/2+1+floor(yi/dx)); 
%-------------------------------------------------------------------------%



%--------------------------Propagagación----------------------------------%
%tamaño del paso de propagacion dz
dz=d/256;
%propagador normalizado
prop=exp(-i*(Kx.^2+Kx'.^2).*dz/4);
% verificando el buen muestreo de la fase del propagador
figure(4)
imagesc(kx,kx,(angle(prop)))
title(['Fase del Propagador'])

% Vector de propagación de 0 a d(hasta distancia de reconstrucción) con paso de dz
zp=0:dz:d;

%calculando el tmaño de zp
sizezp=size(zp);

%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

  %----------------------funcion a propagar----------------------%
 
g=HGo;

gy(:,1)=g(N/2+1,:);
gx(:,1)=g(:,N/2+1);

%ciclo de propagacion

for ii=2:sizezp(2) %corriendo todos los valores de zp a excepcion de z=0
    
    %Transformada de fourier de campo a propagar
    G=fftshift(fft2(g));
    % figure(4)
    % imagesc(x,x,abs(G))
    
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
    
    %cuardando campo transversal
    gy(:,ii)=g(N/2+1,:);
    gx(:,ii)=g(:,N/2+1);
    pxy=g(N,N); g(N,N)=1;
    
    figure(5)
    imagesc(x,x,abs(g))
    title(['Propagación exacta a z=',num2str(zp(ii))])
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    g(N,N)=pxy;
    
    figure(6)
    imagesc(zp,x,abs(gy))
    set(gca,'YDir','normal')
    title('Campo propagado corte en y')
    colormap(mapgreen)
    pause(0.001)
end

%graficando campos tranversales

%cintura de gaussiana
ws=sqrt(1+zp.^2/4);
%cintura de Laguerre
sigmaH=ws*sqrt(2*nu+1);
%Graficando la propagación junto con las normales a los frente de ondas
figure(7)
pcolor(zp,x,abs(gy).^.8)
shading interp
colormap(mapgreen)
title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(zp,-(1/m1)*zp+lo+xt,'y','LineWidth',2.5)
plot(zp,(1/m1)*zp+lo+xt,'c','LineWidth',2.5)
plot(zp,(1/m2)*zp-lo+xt,'y','LineWidth',2.5)
plot(zp,-(1/m2)*zp-lo+xt,'c','LineWidth',2.5)
%graficando la cintura de Laguerre
plot(zp,sigmaH,'b','LineWidth',2.5)
plot(zp,-sigmaH,'b','LineWidth',2.5)
hold off