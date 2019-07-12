%--------------propagacion paraxial por espectro angular------------------%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz

%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Numero de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

Dz=0.75;                    %tamaño de la ventana en z

%indices del Laguerre Gauss
nu=19;mu=0;
%cintura de gaussiana
ws=sqrt(1+Dz.^2/4);
%cintura de Laguerre
sigmaL=ws*sqrt(2*(2*nu+mu+1));

%muestreo de los vectores espaciales
Dx=(2*sigmaL)*1.1;          % Tamano de ventana para el obstáculo centrado en 0
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
%coordenadas polares para el Laguerre
[TH,R]=cart2pol(X,X');


%--------------Laguerre Gauss normalizado(sin unidades físicas)-----------&

%--------------------Laguerre Gauss sin ostrucción------------------------%


%Las dos soluciones del Laguerre Gauss
Ln=exp(i*mu*TH).*LaguerreG(nu,abs(mu),R.^2);
Xn=exp(i*mu*TH).*XLaguerreG(45,nu,abs(mu),R.^2);

%graficando Laguerre
figure(1)
pcolor(x,x,abs(Ln).^.75)
title(['Laguerre con n=',num2str(nu),' m=',num2str(mu)])
axis square
shading flat
colormap(mapgreen)

% Graficando cortes transversales x
figure(2)
subplot(2,1,1);
plot(x,abs(Ln(:,N/2+1)).^2)
title(['Laguerre transversal x con n=',num2str(nu),' m=',num2str(mu)])
subplot(2,1,2);
plot(x,abs(Xn(:,N/2+1)).^2,'r')
title(['XLaguerre transversal x con n=',num2str(nu),' m=',num2str(mu)])


%---------------------Laguerre Gauss con obstrucción----------------------%
lo=2;               %tamaño de la obstrucción (radio de la obstruccion)
xt=0; yt=0;         %tralado de la obstruccion en coordenadas cartesianas
[thetao,ro]=cart2pol(X-xt,X');      %aplicando la traslación a coordendas xy
obo=double(ro<=lo);
clear thetao ro
% Laguerre con obstrucción
Lno=Ln.*(1-obo);
%graficando Laguerre con obstrucción
figure(3)
pcolor(x,x,abs(Lno).^.75)
axis square
shading flat
colormap(mapgreen)
title(['Laguerre con n=',num2str(nu),' m=',num2str(mu),' con obstrucción'])

%-------------Calculando ángulo y distancia de reconstrucción-------------$

phaseH=unwrap(angle(Ln(:,N/2+1)-i*Xn(:,N/2+1)));
wavefrontdr=diff(phaseH')./diff(x);
wavefrontdz=2;                                           % Siempre es 2 (k=2)  
m1=-wavefrontdz/wavefrontdr(floor((lo+xt)/dx)+N/2+1);    % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
m2=-wavefrontdz/wavefrontdr(floor((lo-xt)/dx)+N/2+1);    % Pendiente 2             
d1=(lo+xt)*m1;                              % Distancia de reconstrucción debido a la pendiente 1                            
d2=(lo-xt)*m2;                              % Distancia de reconstrucción debido a la pendiente 2 

% d=max(d1,d2);       %distancia de reconstruccion es el maximo de d1 y d2
d=1;                  %casos extremos y con 2 es la distncia de Rayleigh

%--------------------------Propagagación----------------------------------%
%tamaño del paso de propagacion dz
dz=d/2^6;
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
g=Lno;
%  g=(Ln-i*Xn).*(1-obo);
%cuardando campo transversal en z=0
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
    pxy=g(N,N); g(N,N)=.5;

    figure(5)
    imagesc(x,x,abs(g).^.6)
    title(['Propagación exacta a z=',num2str(zp(ii))])
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    g(N,N)=pxy;
    
    figure(6)
    imagesc(zp,x,abs(gy).^.75)
    title('Campo propagado corte en x')
    colormap(mapgreen)
    set(gca,'YDir','normal')
    pause(0.001)
end

%graficando campos tranversales


%Graficando la propagación junto con las normales a los frente de ondas
figure(7)
pcolor(zp,x,abs(gy).^.6)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en x')
xlabel('z/L_D') % x-axis label
ylabel('x/w(z)') % y-axis label
hold on
%graficando las rectas de las normales al frente de onda
plot(zp,-(1/m1)*zp+lo+xt,'y','LineWidth',2.5)
plot(zp,(1/m1)*zp+lo+xt,'c','LineWidth',2.5)
plot(zp,(1/m2)*zp-lo+xt,'y','LineWidth',2.5)
plot(zp,-(1/m2)*zp-lo+xt,'c','LineWidth',2.5)
% %graficando la cintura de Laguerre
% plot(zp,sigmaL,'b','LineWidth',2.5)
% plot(zp,-sigmaL,'b','LineWidth',2.5)
% plot(zp,ws,'r','LineWidth',2.5)
hold off

figure(8)
pcolor(zp,x,abs(gy).^.6)
shading interp
colormap(mapgreen)
title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(zp,(1/m1)*zp+lo+xt,'c','LineWidth',2.5)
plot(zp,-(1/m2)*zp-lo+xt,'c','LineWidth',2.5)
%graficando la cintura de Laguerre
plot(zp,sigmaL,'b','LineWidth',2.5)
plot(zp,-sigmaL,'b','LineWidth',2.5)
%asintotas de la hiperbola
plot(zp,(sqrt(2*(2*nu+mu+1))/2)*zp,'r','LineWidth',2.5)
plot(zp,-(sqrt(2*(2*nu+mu+1))/2)*zp,'r','LineWidth',2.5)
hold off

figure(9)
pcolor(zp,x,abs(gy).^.6)
shading interp
colormap(mapgreen)
title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(zp,-(1/m1)*zp+lo+xt,'y','LineWidth',2.5)
plot(zp,(1/m2)*zp-lo+xt,'y','LineWidth',2.5)
%graficando la cintura de Laguerre
plot(zp,sigmaL,'b','LineWidth',2.5)
plot(zp,-sigmaL,'b','LineWidth',2.5)
%asintotas de la hiperbola
plot(zp,(sqrt(2*(2*nu+mu+1))/2)*zp,'r','LineWidth',2.5)
plot(zp,-(sqrt(2*(2*nu+mu+1))/2)*zp,'r','LineWidth',2.5)
hold off