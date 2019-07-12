%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Laguerre Gauss------------------------%
nu=19;mu=0;
% Parametros fisicos [micras]
wo=0.02644;
lamb=0.00006328;
k=2*pi/lamb;
zl=k*wo^2/2;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=5*zl; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^6;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=0:dz:Dz;       % Vector z de propagacion

% Cintura de gaussiana a la distancia del tamaño de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2);
% Cintura de Laguerre a la distancia del tamaño de la ventana en z
sigmaL=ws*sqrt((2*nu+mu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tamaño de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(2*sigmaL)*1.3;  % Tamaño de la ventana
dx=Dx/N;            % Resolución
x=n*dx;             % Vector x (plano de campo a propapar)
[X]=meshgrid(x);
% Muestreo del vector de frecuencia
Du=1/dx;        % Tamaño de la ventana de Fourier
du=1/Dx;        % Resolución
u=n*du;         % Vector
[U]=meshgrid(u);
% Vectores kx,ky
kx=2*pi*u;
[Kx]=meshgrid(kx); 
% Coordenadas polares para el Laguerre
[TH,R]=cart2pol(X,X');
%-----------------Laguerre Gauss sin ostrucción en z=0--------------------%

% Las dos soluciones del Laguerre Gauss
Ln=exp(i*mu*TH).*LaguerreG(nu,abs(mu),2*R.^2/wo^2);
Xn=exp(i*mu*TH).*XLaguerreG(45,nu,abs(mu),2*R.^2/wo^2);

% funciones Hankel
H1=Ln+i*Xn;
H2=Ln-i*Xn;

% Función a propagar (Recuerde normalizar respecto al máximo de la función
% en cuestión)
g=Ln;
% Graficando Laguerre
figure(1)
pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^.75)
axis square
shading flat
colormap(mapgreen)
xlabel('x') 
ylabel('y') 

figure(1)
plot(abs(g(N/2+1,:)))

%---------------Funcion a propagar con obstrucción en z=0-----------------%
lo=sqrt(2*nu+mu+1)*wo;%0.057/2;%0.057;               %tamaño de la obstrucción (radio de la obstruccion)
xt=0; yt=0;%0.052; %xt=0.75 %traslado de la obstruccion en coordenadas cartesianas
[thetao,ro]=cart2pol(X-xt,X'-yt);     %aplicando la traslación a coordendas xy
obo=double(ro<=lo);  %creando la obstrucción   
clear thetao ro      %limpienado coordenadas trasladadas
% Función a propagar con obstrucción
% g=g.*(1-obo);
% Recuerde normalizar respecto al máximo de la función en cuestión
pxy=g(1,1); g(1,1)=1.9189; % Para H2 
% Graficando función a propagar con obstrucción
figure(3)
pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)
brighten(.45)
g(1,1)=pxy;

%-------------------------Propagagación física----------------------------%

%propagador paraxial
prop=exp(-i*pi*lamb*dz*(U.^2+(U').^2));
figure(2)
imagesc(u,u,(angle(prop)))
title(['Propagador'])

%calculando el tamaño de zp
sizezp=size(z);

%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

%guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);
gy(:,1)=g(:,N/2+1);


%ciclo de propagacion
% for ii=2:sizezp(2) %corriendo todos los valores de zp a excepcion de z=0
%     
%     %Transformada de fourier de campo a propagar
%     G=fftshift(fft2(g));
%     %obteniendo el campo propagado
%     g=ifft2(fftshift(G.*prop));
%     %cuardando campo transversal
%     gx(:,ii)=g(N/2+1,:);
%     gy(:,ii)=g(:,N/2+1);
% 
%     %campo g propagado
%     fig4=figure(4);
%     fig4.Name=([' z = ',num2str(z(ii))]);
%     pxy=g(N,N); g(N,N)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
%     imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
%     axis square
%     colormap(mapgreen)
%     set(gca,'YDir','normal')
%     g(N,N)=pxy;
%     brighten(.45)
%     axis1=gca;
%     set(axis1,'FontSize',13);
%     xlabel('$x$','Interpreter','latex','FontSize',18) 
%     ylabel('$y$','Interpreter','latex','FontSize',18)
% 
%     %campo g en y,z
%     fig5=figure(5);
% %   Para HL2 (para que estén normalizados HL1 y Hl2)
%     pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
%     imagesc(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
%     fig5.Name=([' z = ',num2str(z(ii))]);
%     colormap(mapgreen)
%     set(gca,'YDir','normal')
%     gy(1,1)=pxy;
%     brighten(0.45)
%     pbaspect([2.5 1 2])
%     axis1=gca;
%     set(axis1,'FontSize',13);
%     xlabel('$z$','Interpreter','latex','FontSize',18) 
%     ylabel('$x$','Interpreter','latex','FontSize',18) 
%     
%     g2=LaguerreG(nu,abs(mu),2*x.^2/(wo^2*(1+(ii*dz)^2/zl^2)));
%     figure(1)
%     plot(abs(g(:,N/2+1)))
%     hold on
%     plot(abs(g2),'r')
%     hold off
%     
%     pause(.01)
% end
% 
% fig9=figure(9);
% %   Para HL2 (para que estén normalizados HL1 y Hl2)
% pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
% pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
% shading interp
% colormap(mapgreen)
% fig9.Name=([' z = ',num2str(z(ii))]);
% gy(1,1)=pxy;
% % title('Campo propagado corte en x')
% axis1=gca;
% set(axis1,'FontSize',13);
% xlabel('$z$','Interpreter','latex','FontSize',18) % x-axis label
% ylabel('$x$','Interpreter','latex','FontSize',18) % y-axis label
% brighten(0.45)
% pbaspect([2.5 1 2])

%-------------Calculando ángulo y distancia de reconstrucción-------------$

% phaseH=unwrap(angle(H2(:,N/2+1)));
% wavefrontdr=diff(phaseH')./diff(x);
% wavefrontdz=k;                                           % derivada de kz 
% m1=wavefrontdz/wavefrontdr(floor((yt+lo)/dx)+N/2+1);    % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
% m2=wavefrontdz/wavefrontdr(floor((yt-lo)/dx)+N/2+1);    % Pendiente 2             
% d1=(lo+yt)*m1;                              % Distancia de reconstrucción debido a la pendiente 1                            
% d2=(lo-yt)*m2;                              % Distancia de reconstrucción debido a la pendiente 2 

% d=max(d1,d2);       %distancia de reconstruccion es el maximo de d1 y d2

    r1=yt+lo;
    r2=yt-lo;
    r3=yt+lo;
    r4=yt-lo;
    
    r11=(yt+lo)/(3/4);
    r13=(yt+lo)/(3/4);

    r21=-(yt+lo)/(3/4);
    r23=-(yt+lo)/(3/4);
    
    Rzi=z(1)+zl^2/z(1);
    wzi=wo*sqrt(1+(z(1)/zl)^2);
    phizi=(2*nu+mu+1)*atan(z(1)/zl);
    % Terminos de haz para todo z
    Rz=z+zl^2./z;
    wz=wo*sqrt(1+(z./zl).^2);
    phiz=(2*nu+mu+1)*atan(z./zl);
    % funcion Hankel 1 en el plano y,z o x,z en z(ii) de propagacion
    H1r=(LaguerreG(nu,abs(mu),2*x.^2./(wzi^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*x.^2./(wzi^2)))...
       .*exp(i*(k*x.^2./(2*Rzi)-phizi));
    % funcion Hankel 1 en el plano y,z o x,z en z(ii) de propagacion
    H2r=(LaguerreG(nu,abs(mu),2*x.^2/(wzi^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*x.^2/(wzi^2)))...
       .*exp(i*(k*x.^2/(2*Rzi)-phizi));
   
    Hz1=(LaguerreG(nu,abs(mu),2*r1.^2./(wz.^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*r1.^2./(wz.^2)))...
       .*exp(i*(k*r1^2./(2*Rz)-phiz));
   
    Hz2=(LaguerreG(nu,abs(mu),2*r2.^2./(wz.^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*r2.^2./(wz.^2)))...
       .*exp(i*(k*r2^2./(2*Rz)-phiz));
    
   
    Hz3=(LaguerreG(nu,abs(mu),2*r3.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r3.^2./(wz.^2)))...
       .*exp(i*(k*r3^2./(2*Rz)-phiz));
   
    Hz4=(LaguerreG(nu,abs(mu),2*r4.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r4.^2./(wz.^2)))...
       .*exp(i*(k*r4^2./(2*Rz)-phiz));
    
    %desdoblando fases
    phaseHr1=unwrap(angle(H1r));    phaseHr2=unwrap(angle(H2r));
    phaseHz1=unwrap(angle(Hz1));    phaseHz2=unwrap(angle(Hz2));
    phaseHz3=unwrap(angle(Hz3));    phaseHz4=unwrap(angle(Hz4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wavefrontdr1=gradient(phaseHr1)/dx;
    wavefrontdr2=gradient(phaseHr2)/dx;
    % Se agrega el k por la derivada del término exp(i*(mu*phi+k*z))
    wavefrontdz1=(gradient(phaseHz1)/dz)+k;
    wavefrontdz2=(gradient(phaseHz2)/dz)+k;
    wavefrontdz3=(gradient(phaseHz3)/dz)+k;
    wavefrontdz4=(gradient(phaseHz4)/dz)+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m1=wavefrontdz1(floor(z(1)/dz+1))./wavefrontdr2(N/2+1+floor(r1/dx));   % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
    m2=wavefrontdz2(floor(z(1)/dz+1))./wavefrontdr2(N/2+1+floor(r2/dx));   % Pendiente 2
    m3=wavefrontdz3(floor(z(1)/dz+1))./wavefrontdr1(N/2+1+floor(r3/dx));
    m4=wavefrontdz4(floor(z(1)/dz+1))./wavefrontdr1(N/2+1+floor(r4/dx));

    m11=wavefrontdz1(floor(z(1)/dz+1))./wavefrontdr2(N/2+1+floor(r11/dx));   % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
    m13=wavefrontdz3(floor(z(1)/dz+1))./wavefrontdr1(N/2+1+floor(r13/dx));

    m21=wavefrontdz1(floor(z(1)/dz+1))./wavefrontdr2(N/2+1+floor(r21/dx));   % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
    m23=wavefrontdz3(floor(z(1)/dz+1))./wavefrontdr1(N/2+1+floor(r23/dx));

%Graficando la propagación junto con las normales a los frente de ondas
figure(10)
pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
shading interp
colormap(mapgreen)
gy(1,1)=pxy;
% title('Campo propagado corte en x')
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$z$','Interpreter','latex','FontSize',18) % x-axis label
ylabel('$x$','Interpreter','latex','FontSize',18) % y-axis label
hold on
%graficando las rectas de las normales al frente de onda (se divide a los (1/m1)*z+xt+lo) por (sqrt(2)*wo) para conservar la relación en las unidades normalizadas)
plot(z/(zl/2),((1/m1)*z+yt+lo)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),(-(1/m1)*z+yt+lo)/(sqrt(2)*wo),'r--','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z+yt-lo)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),(-(1/m2)*z+yt-lo)/(sqrt(2)*wo),'r--','LineWidth',1.5)

plot(z/(zl/2),((1/m11)*z+r11)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),((1/m13)*z+r11)/(sqrt(2)*wo),'r--','LineWidth',1.5)

plot(z/(zl/2),((1/m21)*z+r21)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),((1/m23)*z+r21)/(sqrt(2)*wo),'r--','LineWidth',1.5)

hold off
brighten(0.45)
pbaspect([2.5 1 2])

%---Creando trayectoria de la curva de la sombra por segmentos de rayos--%

% Vectores de los rayos, 2 por cada esquina de la obstruccion en el plano
% y,z o x,z
rx1=zeros(1,sizezp(2));
rx2=zeros(1,sizezp(2));
rx3=zeros(1,sizezp(2));
rx4=zeros(1,sizezp(2));

% Primer posición del vector en z=0 para cada uno de las esquinas de la
% sombra justo en la obstrucción
rx1(1)=yt+lo;
rx2(1)=yt-lo;
rx3(1)=yt+lo;
rx4(1)=yt-lo;

% posición de la sombra en la obstrucción en z=0
r1=yt+lo;
r2=yt-lo;
r3=yt+lo;
r4=yt-lo;

% las pendientes ya las calculamos, 2 por cada esquina 
m3=-m1;
m4=-m2;
    
for ii=2:sizezp(2) %corriendo todos los valores de sp
    
    % Dado z encontramos z+dz y la posición de r que le corresponde a este
    % punto z+dz con la pendiente del frente de onda que hay entre estos 
    % dos puntos
    r1=(1/m1)*(z(ii)-z(ii-1))+r1;  
    r2=(1/m2)*(z(ii)-z(ii-1))+r2; 
    r3=(1/m3)*(z(ii)-z(ii-1))+r3;
    r4=(1/m4)*(z(ii)-z(ii-1))+r4;
    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta
    rx1(ii)=r1;
    rx2(ii)=r2;
    rx3(ii)=r3;
    rx4(ii)=r4;
    %-----------Calculando las pendientes del siguiente paso--------------%
    % Terminos de haz a la distancia z(ii)
    Rzi=z(ii)+zl^2/z(ii);
    wzi=wo*sqrt(1+(z(ii)/zl)^2);
    phizi=(2*nu+mu+1)*atan(z(ii)/zl);
    % Terminos de haz para todo z
    Rz=z+zl^2./z;
    wz=wo*sqrt(1+(z./zl).^2);
    phiz=(2*nu+mu+1)*atan(z./zl);
    % funcion Hankel 1 en el plano y,z o x,z en z(ii) de propagacion
    H1r=(LaguerreG(nu,abs(mu),2*x.^2./(wzi^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*x.^2./(wzi^2)))...
       .*exp(i*(k*x.^2./(2*Rzi)-phizi));
    % funcion Hankel 1 en el plano y,z o x,z en z(ii) de propagacion
    H2r=(LaguerreG(nu,abs(mu),2*x.^2/(wzi^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*x.^2/(wzi^2)))...
       .*exp(i*(k*x.^2/(2*Rzi)-phizi));
   
    % El if es para que cuando el rayo cruce el eje el rayo se rija ahora
    % por H1 y no por H2
    if r1>0
    Hz1=(LaguerreG(nu,abs(mu),2*r1.^2./(wz.^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*r1.^2./(wz.^2)))...
       .*exp(i*(k*r1^2./(2*Rz)-phiz));
    else
    Hz1=(LaguerreG(nu,abs(mu),2*r1.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r1.^2./(wz.^2)))...
       .*exp(i*(k*r1^2./(2*Rz)-phiz));
    end
    
    % El if es para que cuando el rayo cruce el eje el rayo se rija ahora
    % por H1 y no por H2
    if r2<0
       Hz2=(LaguerreG(nu,abs(mu),2*r2.^2./(wz.^2))...
       -i*XLaguerreG(45,nu,abs(mu),2*r2.^2./(wz.^2)))...
       .*exp(i*(k*r2^2./(2*Rz)-phiz));
    else
       Hz2=(LaguerreG(nu,abs(mu),2*r2.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r2.^2./(wz.^2)))...
       .*exp(i*(k*r2^2./(2*Rz)-phiz));
    end
    % Les pongo el signo menos porque parece que así se calcula más fácil y
    % al final le cambio el signo a la derivada
    Hz3=(LaguerreG(nu,abs(mu),2*r3.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r3.^2./(wz.^2)))...
       .*exp(i*(k*r3^2./(2*Rz)-phiz));
   
    Hz4=(LaguerreG(nu,abs(mu),2*r4.^2./(wz.^2))...
       +i*XLaguerreG(45,nu,abs(mu),2*r4.^2./(wz.^2)))...
       .*exp(i*(k*r4^2./(2*Rz)-phiz));
    
    %desdoblando fases
    phaseHr1=unwrap(angle(H1r));    phaseHr2=unwrap(angle(H2r));
    phaseHz1=unwrap(angle(Hz1));    phaseHz2=unwrap(angle(Hz2));
    phaseHz3=unwrap(angle(Hz3));    phaseHz4=unwrap(angle(Hz4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wavefrontdr1=gradient(phaseHr1)/dx;
    wavefrontdr2=gradient(phaseHr2)/dx;
    % Se agrega el k por la derivada del término exp(i*(mu*phi+k*z))
    wavefrontdz1=(gradient(phaseHz1)/dz)+k;
    wavefrontdz2=(gradient(phaseHz2)/dz)+k;
    wavefrontdz3=(gradient(phaseHz3)/dz)+k;
    wavefrontdz4=(gradient(phaseHz4)/dz)+k;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Este If es para que cuando el rayo haya cruzado el eje se use la
    % derivada de H1 respecto a r 
    if r1>0
        m1=wavefrontdz1(floor(z(ii)/dz+1))./wavefrontdr2(N/2+1+floor(r1/dx));   % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
    else
        m1=wavefrontdz1(floor(z(ii)/dz+1))./wavefrontdr1(N/2+1+floor(r1/dx));
    end
    if r2<0
        m2=wavefrontdz2(floor(z(ii)/dz+1))./wavefrontdr2(N/2+1+floor(r2/dx));   % Pendiente 2
    else
        m2=wavefrontdz2(floor(z(ii)/dz+1))./wavefrontdr1(N/2+1+floor(r2/dx));
    end
    m3=wavefrontdz3(floor(z(ii)/dz+1))./wavefrontdr1(N/2+1+floor(r3/dx));
    m4=wavefrontdz4(floor(z(ii)/dz+1))./wavefrontdr1(N/2+1+floor(r4/dx));
end
%
figure(11)
pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
pcolor(z/(zl/2),x/(sqrt(2)*wo),angle(gy))
shading interp
% colormap(mapgreen)
gy(1,1)=pxy;
% title('Campo propagado corte en x')
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$z$','Interpreter','latex','FontSize',18) % x-axis label
ylabel('$x$','Interpreter','latex','FontSize',18) % y-axis label
hold on
%graficando las rectas de las normales al frente de onda (se divide a las rxj por (sqrt(2)*wo) para conservar la relación en las unidades normalizadas)
plot(z/(zl/2),rx1/(sqrt(2)*wo),'c','LineWidth',1.5) % 'c--'
plot(z/(zl/2),rx2/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),rx3/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),rx4/(sqrt(2)*wo),'r','LineWidth',1.5)
hold off
% brighten(0.45)
pbaspect([2.5 1 2])

% Quiver
r=3.12:0.005:3.14;
z00=0.29:0.01:0.31;z0=z00';
ddr=diff(r);ddz=diff(z00);

Rzi=z(ii)+zl^2/z(ii);
wzi=wo*sqrt(1+(z(ii)/zl)^2);
phizi=(2*nu+mu+1)*atan(z(ii)/zl);
% Terminos de haz para todo z
Rz=z+zl^2./z;
wz=wo*sqrt(1+(z./zl).^2);
phiz=(2*nu+mu+1)*atan(z./zl);
phangle=2*z0+unwrap(angle((LaguerreG(19,0,2*r.^2./(2*(1+(z0/2).^2)))...
    +i*XLaguerreG(45,19,0,2*r.^2./(2*(1+(z0/2).^2))))));
% phatan=2*z0+atan(XLaguerreG(45,19,0,2*r.^2./(2*(1+(z0/2).^2)))...
%     ./LaguerreG(19,0,2*r.^2./(2*(1+(z0/2).^2))));
[px1,py1]=gradient(phangle);
% [px2,py2]=gradient(phatan);
figure(13)
contour(r,z0,phangle)
axis square
hold on
quiver(r,z0,px1/ddr(1),py1/ddz(1))
hold off
% figure(14)
% contour(r,z0,phatan)
% axis square
% hold on
% quiver(r,z0,px2/ddr(1),py2/ddz(1))
% hold off