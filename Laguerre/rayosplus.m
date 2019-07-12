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
Dz=zl/2; % Tamaño de la ventana en z (distancia a la cual propagar)
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
Ln=exp(i*mu*TH).*laguerregz(nu,mu,wo,zl,R,0); 
% Ln=exp(i*mu*TH).*LaguerreG(nu,abs(mu),2*R.^2/wo^2);
Xn=exp(i*mu*TH).*xlaguerregz(nu,mu,wo,zl,R,0);

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

figure(2)
plot(abs(g(N/2+1,:)))

%---------------Funcion a propagar con obstrucción en z=0-----------------%
sigmaLo=wo*sqrt((2*nu+mu+1));

lo=sigmaLo/4;%0.057/2;%0.057;               %tamaño de la obstrucción (radio de la obstruccion)
xt=0; yt=0;%0.052; %xt=0.75 %traslado de la obstruccion en coordenadas cartesianas
[thetao,ro]=cart2pol(X-xt,X'-yt);     %aplicando la traslación a coordendas xy
obo=double(ro<=lo);  %creando la obstrucción   
clear thetao ro      %limpiando coordenadas trasladadas
% Función a propagar con obstrucción
g=g.*(1-obo);
% Recuerde normalizar respecto al máximo de la función en cuestión
% pxy=g(1,1); g(1,1)=1.9189; % Para H2 
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

%-------------------------Propagagación física----------------------------%
%propagador paraxial
prop=exp(-i*pi*lamb*dz*(U.^2+(U').^2));
figure(4)
imagesc(u,u,(angle(prop)))
title(['Propagador'])
%calculando el tamaño de zp
sizezp=size(z);
%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));
%guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);
gy(:,1)=g(:,N/2+1);
% ciclo de propagacion
for ii=2:size(z,2) %corriendo todos los valores de zp a excepcion de z=0
    
    G=fftshift(fft2(g));
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
%    guardando campo analitico de Laguerre para todo z
%     g=laguerregz(nu,mu,wo,zl,R,z(ii));
%     cuardando campo transversal
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);

    %campo g propagado
    fig5=figure(5);
    fig5.Name=([' z = ',num2str(z(ii))]);
    pxy=g(N,N); g(N,N)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    g(N,N)=pxy;
    brighten(.45)
    axis1=gca;
    set(axis1,'FontSize',13);
    xlabel('$x$','Interpreter','latex','FontSize',18) 
    ylabel('$y$','Interpreter','latex','FontSize',18)

    %campo g en y,z
    fig6=figure(6);
%   Para HL2 (para que estén normalizados HL1 y Hl2)
    pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
    imagesc(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
    fig6.Name=([' z = ',num2str(z(ii))]);
    colormap(mapgreen)
    set(gca,'YDir','normal')
    gy(1,1)=pxy;
    brighten(0.45)
    pbaspect([2.5 1 2])
    axis1=gca;
    set(axis1,'FontSize',13);
    xlabel('$z$','Interpreter','latex','FontSize',18) 
    ylabel('$x$','Interpreter','latex','FontSize',18) 

    pause(.01)
end
figure(7)
imagesc(phase_unwrap(angle(gy)))

fig9=figure(9);
%   Para HL2 (para que estén normalizados HL1 y Hl2)
pxy=gy(1,1); gy(1,1)=1.9189; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
shading interp
colormap(mapgreen)
fig9.Name=([' z = ',num2str(z(ii))]);
gy(1,1)=pxy;
% title('Campo propagado corte en x')
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$z$','Interpreter','latex','FontSize',18) % x-axis label
ylabel('$x$','Interpreter','latex','FontSize',18) % y-axis label
brighten(0.45)
pbaspect([2.5 1 2])
%L=sin(3*pi/8)/sin(pi/4)r
%-------------Calculando ángulo y distancia de reconstrucción-------------$
%posición del rayo en la coordenada radial para los ejes en la cintura
rx=lo;
ry=lo;

%calculo de las pendientes en el punto (sigmaL,0,z)
% para z=cte;
H1r=laguerregz(nu,mu,wo,zl,x,z(1))+1i*xlaguerregz(nu,mu,wo,zl,x,z(1));
% para r=cte;
H1z=laguerregz(nu,mu,wo,zl,rx,z)+1i*xlaguerregz(nu,mu,wo,zl,rx,z);

fr=unwrap(angle(H1r));
fz=unwrap(angle(H1z));

[mzro]=gradientrz(fr,fz,k,dx,dz,rx,0); 
 
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
plot(z/(zl/2),((1/mzro)*z+yt+lo)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),(-(1/mzro)*z+yt+lo)/(sqrt(2)*wo),'r--','LineWidth',1.5)
plot(z/(zl/2),((1/mzro)*z+yt-lo)/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),(-(1/mzro)*z+yt-lo)/(sqrt(2)*wo),'r--','LineWidth',1.5)
hold off
pbaspect([2.5 1 2])

%---Creando trayectoria de la curva de la sombra por segmentos de rayos--%

% Vectores de los rayos, 2 por cada esquina de la obstruccion en el plano
% y,z o x,z
rx1=zeros(1,sizezp(2));
rx2=zeros(1,sizezp(2));

% Primer posición del vector en z=0 para cada uno de las esquinas de la
% sombra justo en la obstrucción
rx1(1)=rx;
rx2(1)=rx;
r1=rx;
r2=rx;
mzr1=mzro;
mzr2=-mzro;

for ii=2:sizezp(2) %corriendo todos los valores de sp
    
    % Dado z encontramos z+dz y la posición de r que le corresponde a este
    % punto z+dz con la pendiente del frente de onda que hay entre estos 
    % dos puntos
    r1=(1/mzr1)*(z(ii)-z(ii-1))+r1;
    r2=(1/mzr2)*(z(ii)-z(ii-1))+r2;  

    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta
    rx1(ii)=r1;
    rx2(ii)=r2;
    %-----------Calculando las pendientes del siguiente paso--------------%
    %calculo de la pendiente en (rx,z(ii)) de Hankel 1
    % para z=cte;
    H1r1=laguerregz(nu,mu,wo,zl,x,z(ii))+1i*xlaguerregz(nu,mu,wo,zl,x,z(ii));
    % para r=cte;
    H1z1=laguerregz(nu,mu,wo,zl,r1,z)+1i*xlaguerregz(nu,mu,wo,zl,r1,z);
    fr=unwrap(angle(H1r1));
    fz=unwrap(angle(H1z1));
    [mzr1]=gradientrz(fr,fz,k,dx,dz,r1,z(ii)); 
    
    %calculo de la pendiente en (rx,z(ii)) de Hankel 2
    %Aqui el rayo que entra debido a Hankel 1 despues se rige por Hankel 2
    %cuando pasa por el origen
    if r2<0
    % para z=cte;
    H1r2=laguerregz(nu,mu,wo,zl,x,z(ii))+1i*xlaguerregz(nu,mu,wo,zl,x,z(ii));
    % para r=cte;
    H1z2=laguerregz(nu,mu,wo,zl,r2,z)+1i*xlaguerregz(nu,mu,wo,zl,r2,z);
    fr=unwrap(angle(H1r2));
    fz=unwrap(angle(H1z2));
    [mzr2]=gradientrz(fr,fz,k,dx,dz,r2,z(ii)); 
    else
    % para z=cte;
    H2r2=laguerregz(nu,mu,wo,zl,x,z(ii))-1i*xlaguerregz(nu,mu,wo,zl,x,z(ii));
    % para r=cte;
    H2z2=laguerregz(nu,mu,wo,zl,r2,z)-1i*xlaguerregz(nu,mu,wo,zl,r2,z);
    fr=unwrap(angle(H2r2));
    fz=unwrap(angle(H2z2));
    [mzr2]=gradientrz(fr,fz,k,dx,dz,r2,z(ii)); 
    end

end
%
figure(11)
imagesc(z/(zl/2),x/(sqrt(2)*wo),abs((gy)))
% colormap(mapgreen)
% title('Campo propagado corte en x')
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$z$','Interpreter','latex','FontSize',18) % x-axis label
ylabel('$x$','Interpreter','latex','FontSize',18) % y-axis label
hold on
%graficando las rectas de las normales al frente de onda (se divide a las rxj por (sqrt(2)*wo) para conservar la relación en las unidades normalizadas)
plot(z/(zl/2),rx1/(sqrt(2)*wo),'r','LineWidth',1.5) % 'c--'
plot(z/(zl/2),rx2/(sqrt(2)*wo),'c','LineWidth',1.5) % 'c--'
plot(z/(zl/2),-rx1/(sqrt(2)*wo),'r','LineWidth',1.5) % 'c--'
plot(z/(zl/2),-rx2/(sqrt(2)*wo),'c','LineWidth',1.5) % 'c--'
hold off
% brighten(0.45)
pbaspect([2.5 1 2])

