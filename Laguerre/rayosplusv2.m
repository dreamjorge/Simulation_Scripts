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
Dz=zl; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^10;          % Número de puntos en z
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
Ln=exp(1i*mu*TH).*laguerregz(nu,mu,wo,zl,R,0); 
% Ln=exp(i*mu*TH).*LaguerreG(nu,abs(mu),2*R.^2/wo^2);
Xn=exp(1i*mu*TH).*xlaguerregz(nu,mu,wo,zl,R,0);

% funciones Hankel
H1=Ln+1i*Xn;
H2=Ln-1i*Xn;

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
pxy=max(max(g));

figure(2)
plot(abs(g(N/2+1,:)))

%---------------Funcion a propagar con obstrucción en z=0-----------------%
sigmaLo=wo*sqrt((2*nu+mu+1));

lo=sigmaLo/4;       %tamaño de la obstrucción (radio de la obstruccion)
[~,ro]=cart2pol(X,X');     %aplicando la traslación a coordendas xy
obo=double(ro<=lo);  %creando la obstrucción   
clear thetao ro      %limpiando coordenadas trasladadas
% Función a propagar con obstrucción
g=g.*(1-obo);

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

%---------------calculo de los rayos en el punto (rx,z=0)-----------------%
%posición del rayo en la coordenada radial para los ejes en la cintura
rx=lo; %en la parte radial de la obstruccion
%calculo de las pendientes en el punto (rx,z) para Hankel 1
% para z=cte;
H1r=laguerregz(nu,mu,wo,zl,x,z(1))+1i*xlaguerregz(nu,mu,wo,zl,x,z(1));
% para r=cte;
H1z=laguerregz(nu,mu,wo,zl,rx,z)+1i*xlaguerregz(nu,mu,wo,zl,rx,z);
% fases
fr=unwrap(angle(H1r));
fz=unwrap(angle(H1z));
% pendiente
[mzro]=gradientrz(fr,fz,k,dx,dz,rx,0); 
%---Creando trayectoria de la curva de la sombra por segmentos de rayos--%
% Vectores de los rayos, 2 por cada esquina de la obstruccion(en la parte radial)
rx1=zeros(1,size(z,2));
rx2=zeros(1,size(z,2));
% Primer posición del vector en z=0 para cada uno de las esquinas de la
% sombra justo en la obstrucción
rx1(1)=rx;
rx2(1)=rx;
r1=rx;
r2=rx;
mzr1=mzro;
mzr2=-mzro;
%-----------------------------fin de rayos--------------------------------%
%-------------------------Propagagación física----------------------------%
%propagador paraxial
prop=exp(-1i*pi*lamb*dz*(U.^2+(U').^2));
figure(4)
imagesc(u,u,(angle(prop)))
title('Propagador')
%calculando el tamaño de zp
sizezp=size(z);
%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));
%guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);
gy(:,1)=g(:,N/2+1);


% vidObj1 = VideoWriter('LGo.avi');
% vidObj1.Quality = 100;
% vidObj1.FrameRate =30;
% open(vidObj1);

or = [1 .5 0];
for ii=2:sizezp(2) %corriendo todos los valores de sp
 %---------------------------Calculo de Rayos----------------------------%   
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
 %-----------------------Fin de Calculo de Rayos--------------------------%   
   
    %propagacion del campo
    G=fftshift(fft2(g));
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
%   cuardando campo transversal
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);
    %campo g propagado
    pxyz=g(1,1);
    g(1,1)=pxy;
    figure(5);
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^.75)
    hold on
    axis square
    colormap(mapgreen )
    set(gca,'YDir','normal')
    axis1=gca;
    set(axis1,'FontSize',13);
    xlabel('$x$','Interpreter','latex','FontSize',18) 
    ylabel('$y$','Interpreter','latex','FontSize',18)
    plot(r1/(sqrt(2)*wo),0,'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(0,r1/(sqrt(2)*wo),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(r1/(sqrt(2)*wo)*cos(pi/4),r1/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(r2/(sqrt(2)*wo),0,'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(0,r2/(sqrt(2)*wo),'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(r2/(sqrt(2)*wo)*cos(pi/4),r2/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(-r1/(sqrt(2)*wo),0,'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(0,-r1/(sqrt(2)*wo),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(-r1/(sqrt(2)*wo)*cos(pi/4),-r1/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(-r2/(sqrt(2)*wo),0,'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(0,-r2/(sqrt(2)*wo),'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(-r2/(sqrt(2)*wo)*cos(pi/4),-r2/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(r1/(sqrt(2)*wo)*cos(pi/4),-r1/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(r2/(sqrt(2)*wo)*cos(pi/4),-r2/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','y')
    plot(-r1/(sqrt(2)*wo)*cos(pi/4),r1/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','r')
    plot(-r2/(sqrt(2)*wo)*cos(pi/4),r2/(sqrt(2)*wo)*sin(pi/4),'+','MarkerSize',10,'LineWidth',2,'color','y')
    text(-7.5,7.5,[' z = ',num2str(z(ii)/zl)],'Color','red','FontSize',12)
    hold off
    g(1,1)=pxyz;
    pause(.1)
%     writeVideo(vidObj1, getframe(gca));

end
% close(vidObj1);


