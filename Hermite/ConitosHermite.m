%--------------propagacion paraxial por espectro angular------------------%
mapgreen = AdvancedColormap('kg',256,[0 255]/255);  %color del haz

%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Numero de puntos
n=-N/2:N/2-1;               % Vector indicial igualmente espaciado

%muestreo de los vectores espaciales
wo=.5;
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
[HGx,NHGx]= hermitegauss(nu,x/wo);
[HGy,NHGy]= hermitegauss(mu,x/wo);

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
imagesc(x,x,angle(HG))
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


%---------------------Laguerre Gauss con obstrucción----------------------%
lo=2;               %tamaño de la obstrucción (radio de la obstruccion)
xt=2; yt=0;         %tralado de la obstruccion en coordenadas cartesianas
ob=(abs(X'-yt)<=(lo)).*(abs(X-xt)<=(lo));


% Hermite con obstrucción
HGo=HG.*(1-ob);
%graficando Hermite con obstrucción
figure(3)
imagesc(x,x,abs(HGo))
axis square
colormap(mapgreen)
title(['Hermite Gauss con n=',num2str(nu),' m=',num2str(mu),' con obstrucción'])


%-------------Calculando ángulo y distancia de reconstrucción-------------$
H=(H1y')*(H1x);


g=angle(H(:,N/2+1));
phasedr=diff(g')./diff(x);
phasedz=2;       % Siempre es 2 (k=2)  
m1=phasedz/phasedr(floor((lo+xt)/dx)+N/2+1);    % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
m2=phasedz/phasedr(floor((lo-xt)/dx)+N/2+1);    % Pendiente 2             
d1=(lo+xt)*m1;                              % Distancia de reconstrucción debido a la pendiente 1                            
d2=(lo-xt)*m2;                              % Distancia de reconstrucción debido a la pendiente 2 

d=max(d1,d2);       %distancia de reconstruccion es el maximo de d1 y d2
d=1;               %casos extremos y con 2 es la distncia de Rayleigh

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