%Hermite analiticos
clear all
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite-----------------------%
nu=20;mu=20;
% Parametros fisicos [micras]
wo=300;
lamb=0.6328;
k=2*pi/lamb;
zl=k*wo^2/2;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2:N/2-1;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl/2;%/5;   % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^6;        % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=(0:Nz-1)*dz;   % Vector z de propagacion
% calculando el tamaño de zp
sizezp=size(z);
% Cintura de gaussiana a la distancia del tamaño de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2)/sqrt(2);
% Cintura de Laguerre a la distancia del tamaño de la ventana en z
niu=max(nu,mu);
sigmaH=ws*sqrt((2*niu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tamaño de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(sigmaH)*2.5;%*1.6/(1.97)^2;  % Tamaño de la ventana
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
%-------------------Cintura de Hermite para todo z------------------------%
wz=wo*sqrt(1+z.^2/zl^2);
sigmaHx=sqrt(2*nu+1)*wz/sqrt(2);
sigmaHy=sqrt(2*mu+1)*wz/sqrt(2);
sigmaH=sigmaHx.^2+sigmaHy.^2;
%------------------------campo inicial en z=0-----------------------------%
[HGx,NHGx]=hermiteg(nu,sqrt(2)*x/wo);
[HGy,NHGy]=hermiteg(mu,sqrt(2)*x/wo);
HG=(HGy')*(HGx);
figure(1)
imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(HG).^2)
axis square
% axis([-3,3,-3,3])
colormap(mapgreen)
set(gca,'YDir','normal')
axis1=gca;
set(axis1,'FontSize',28);
xlabel('$x$','Interpreter','latex','FontSize',28) 
ylabel('$y$','Interpreter','latex','FontSize',28)
hold on
% Dibujando cuadro de cintura de Hermite
rectangle('Position',[-sigmaHx(1)/(sqrt(2)*wo) -sigmaHy(1)/(sqrt(2)*wo) 2*sigmaHx(1)/(sqrt(2)*wo) 2*sigmaHy(1)/(sqrt(2)*wo)],'EdgeColor','b','LineStyle','--','LineWidth',2)
hold off

%obstruccion
lx=sigmaHx(1)/4;
ly=sigmaHy(1)/4;
obx=double(abs(x)<=lx);
oby=double(abs(x)<=ly);
ob=(oby')*obx;
HGo=HG.*(1-ob);

figure(2)
imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(HGo).^2)
axis square
% axis([-3,3,-3,3])
colormap(mapgreen)
set(gca,'YDir','normal')
axis1=gca;
set(axis1,'FontSize',28);
xlabel('$x$','Interpreter','latex','FontSize',28) 
ylabel('$y$','Interpreter','latex','FontSize',28)
%calculamos las pendientes en las esquinas de la obstrucción--------------%
% calculando la pendiente de H11 en  x=lx,y=ly,z=0
xi=lx; yi=ly; zi=0;
HH11yz=hankelH2(1,1,nu,mu,wo,zl,x,yi,zi);
HH11xz=hankelH2(1,1,nu,mu,wo,zl,xi,x,zi);
HH11xy=hankelH2(1,1,nu,mu,wo,zl,xi,yi,z);

%componentes de la fase
fyz=unwrap(angle(HH11yz));
fxz=unwrap(angle(HH11xz));
fxy=unwrap(angle(HH11xy));

%calculo de las pendientes
[mzx11o,mzy11o,mxy11o]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 

% calculando de la pendiente de H12 en  x=lx,y=ly,z=0
HH12yz=hankelH2(1,2,nu,mu,wo,zl,x,yi,zi);
HH12xz=hankelH2(1,2,nu,mu,wo,zl,xi,x,zi);
HH12xy=hankelH2(1,2,nu,mu,wo,zl,xi,yi,z);

%componentes de la fase
fyz=unwrap(angle(HH12yz));
fxz=unwrap(angle(HH12xz));
fxy=unwrap(angle(HH12xy));

%calculo de las pendientes
[mzx12o,mzy12o,mxy12o]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 

% calculando de la pendiente de H21 en  x=lx,y=ly,z=0
HH21yz=hankelH2(2,1,nu,mu,wo,zl,x,yi,zi);
HH21xz=hankelH2(2,1,nu,mu,wo,zl,xi,x,zi);
HH21xy=hankelH2(2,1,nu,mu,wo,zl,xi,yi,z);

%componentes de la fase
fyz=unwrap(angle(HH21yz));
fxz=unwrap(angle(HH21xz));
fxy=unwrap(angle(HH21xy));

%calculo de las pendientes
[mzx21o,mzy21o,mxy21o]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 

% calculando el gradiente de H22 en  x=lx,y=ly,z=0
HH22yz=hankelH2(2,2,nu,mu,wo,zl,x,yi,zi);
HH22xz=hankelH2(2,2,nu,mu,wo,zl,xi,x,zi);
HH22xy=hankelH2(2,2,nu,mu,wo,zl,xi,yi,z);

%componentes de la fase
fyz=unwrap(angle(HH22yz));
fxz=unwrap(angle(HH22xz));
fxy=unwrap(angle(HH22xy));

%calculo de las pendientes
[mzx22o,mzy22o,mxy22o]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 

% Vectores para los rayos en x e y para cada una de las trayectorias de las
% Hankel
% rayos para x
rxz11=zeros(1,sizezp(2)); rxz12=zeros(1,sizezp(2));
rxz21=zeros(1,sizezp(2)); rxz22=zeros(1,sizezp(2));
% rayos para y
ryz11=zeros(1,sizezp(2)); ryz12=zeros(1,sizezp(2));   
ryz21=zeros(1,sizezp(2)); ryz22=zeros(1,sizezp(2));   

% primer posicion del rayo en x
rx11=lx; rx12=lx; rx21=lx; rx22=lx;

% primer posicion del rayo en y
ry11=ly; ry12=ly; ry21=ly; ry22=ly;

% pendientes en z=0 en el punto(lx,ly)
% pendiente de z con respecto a x
mzx11=mzx11o; mzx12=mzx12o; mzx21=mzx21o; mzx22=mzx22o;
% pendiente de z con respecto a y
mzy11=mzy11o; mzy12=mzy12o; mzy21=mzy21o; mzy22=mzy22o;

% Primer valor de los rayos para x
rxz11(1)=rx11; rxz12(1)=rx12; rxz21(1)=rx21; rxz22(1)=rx22;
% Primer valor de los
ryz11(1)=ry11; ryz12(1)=ry12; ryz21(1)=ry21; ryz22(1)=ry22;   


%propagador
prop=exp(-1i*pi*lamb*dz*(U.^2+(U').^2));
figure(3)
imagesc((angle(prop)))
%funcion a propagar
g=HGo;

%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

gy(:,1)=g(N/2+1,:);
gx(:,1)=g(:,N/2+1);

for ii=2:sizezp(2) %corriendo todos los valores de zp
    
    %Transformada de fourier de campo a propagar
    G=fftshift(fft2(g));
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
    gx(:,ii)=g(:,N/2+1);
    gy(:,ii)=g(N/2+1,:);

    % para x
    rx11=(1/mzx11)*(z(ii)-z(ii-1))+rx11;
    rx12=(1/mzx12)*(z(ii)-z(ii-1))+rx12;
    rx21=(1/mzx21)*(z(ii)-z(ii-1))+rx21;
    rx22=(1/mzx22)*(z(ii)-z(ii-1))+rx22;
    % para y
    ry11=(1/mzy11)*(z(ii)-z(ii-1))+ry11;
    ry12=(1/mzy12)*(z(ii)-z(ii-1))+ry12;
    ry21=(1/mzy21)*(z(ii)-z(ii-1))+ry21;
    ry22=(1/mzy22)*(z(ii)-z(ii-1))+ry22;
    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta para cada una de las Hankel
    % para x
    rxz11(ii)=rx11; rxz12(ii)=rx12; rxz21(ii)=rx21; rxz22(ii)=rx22;
    % para y
    ryz11(ii)=ry11; ryz12(ii)=ry12; ryz21(ii)=ry21; ryz22(ii)=ry22;

    %---------------Pendiente de H11 en x=rx11,y=ry11,z=z(ii)-------------%
    xi=rx11; yi=ry11; zi=z(ii);
    HH11yz=hankelH2(1,1,nu,mu,wo,zl,x,yi,zi);
    HH11xz=hankelH2(1,1,nu,mu,wo,zl,xi,x,zi);
    HH11xy=hankelH2(1,1,nu,mu,wo,zl,xi,yi,z);

    %componentes de la fase
    fyz=unwrap(angle(HH11yz));
    fxz=unwrap(angle(HH11xz));
    fxy=unwrap(angle(HH11xy));

    %calculo de las pendientes
    [mzx11,mzy11,mxy11]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 
    %---------------Pendiente de H12 en x=rx12,y=ry12,z=z(ii)-------------%
    xi=rx12; yi=ry12; zi=z(ii);
    HH12yz=hankelH2(1,2,nu,mu,wo,zl,x,xi,zi);
    HH12xz=hankelH2(1,2,nu,mu,wo,zl,xi,x,zi);
    HH12xy=hankelH2(1,2,nu,mu,wo,zl,xi,xi,z);

    %componentes de la fase
    fyz=unwrap(angle(HH12yz));
    fxz=unwrap(angle(HH12xz));
    fxy=unwrap(angle(HH12xy));

    %calculo de las pendientes
    [mzx12,mzy12,mxy12]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 
    %---------------Pendiente de HH21 en x=rx21,y=rx21,z=z(ii)------------%
    xi=rx21; yi=ry21; zi=z(ii);
    HH21yz=hankelH2(2,1,nu,mu,wo,zl,x,yi,zi);
    HH21xz=hankelH2(2,1,nu,mu,wo,zl,xi,x,zi);
    HH21xy=hankelH2(2,1,nu,mu,wo,zl,xi,yi,z);

    %componentes de la fase
    fyz=unwrap(angle(HH21yz));
    fxz=unwrap(angle(HH21xz));
    fxy=unwrap(angle(HH21xy));

    %calculo de las pendientes
    [mzx21,mzy21,mxy21]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 
    %---------------Pendiente de HH22 en x=rx22,y=rx22,z=z(ii)------------%
    xi=rx22; yi=ry22; zi=z(ii);
    HH22yz=hankelH2(2,2,nu,mu,wo,zl,x,yi,zi);
    HH22xz=hankelH2(2,2,nu,mu,wo,zl,xi,x,zi);
    HH22xy=hankelH2(2,2,nu,mu,wo,zl,xi,yi,z);

    %componentes de la fase
    fyz=unwrap(angle(HH22yz));
    fxz=unwrap(angle(HH22xz));
    fxy=unwrap(angle(HH22xy));

    %calculo de las pendientes
    [mzx22,mzy22,mxy22]=gradientxyz(fyz,fxz,fxy,k,dx,dx,dz,xi,yi,zi); 

    figure(10)
    fig10=figure(10);
    fig10.Name=([' z = ',num2str(z(ii))]);
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis1=gca;
    set(axis1,'FontSize',28);
    xlabel('$x$','Interpreter','latex','FontSize',28) 
    ylabel('$y$','Interpreter','latex','FontSize',28)
    hold on
    %cuadro de la cintura en del Hermite en z
    rectangle('Position',[-sigmaHx(ii)/(sqrt(2)*wo) -sigmaHy(ii)/(sqrt(2)*wo) 2*sigmaHx(ii)/(sqrt(2)*wo) 2*sigmaHy(ii)/(sqrt(2)*wo)],'EdgeColor','b','LineStyle','--','LineWidth',2)
    %Propagacion de las cuatro hankel en la esquina superior derecha
    plot(rx11/(sqrt(2)*wo),ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(rx12/(sqrt(2)*wo),ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx21/(sqrt(2)*wo),ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx22/(sqrt(2)*wo),ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %Propagacion de las cuatro Hankel en la esquina superior izquierda
    plot(-rx11/(sqrt(2)*wo),ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(-rx12/(sqrt(2)*wo),ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx21/(sqrt(2)*wo),ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx22/(sqrt(2)*wo),ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %Propagacion de las cuatro Hankel en la esquina inferior izquierna
    plot(-rx11/(sqrt(2)*wo),-ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(-rx12/(sqrt(2)*wo),-ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx21/(sqrt(2)*wo),-ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx22/(sqrt(2)*wo),-ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %Propagacion de las cuatro Hankel en la esquina inferior derecha
    plot(rx11/(sqrt(2)*wo),-ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(rx12/(sqrt(2)*wo),-ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx21/(sqrt(2)*wo),-ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx22/(sqrt(2)*wo),-ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    hold off
    pause(.001)
%     writeVideo(vidObj1, getframe(gca));
end  
% close(vidObj1);

figure(8)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
shading interp
colormap(mapgreen)
% title('Campo propagado corte en y')
hold on
%graficando las rectas de las normales al frente de onda
%graficando los rayos de las normales al frente de onda
plot(z/(zl/2),(rxz11)/(sqrt(2)*wo),'r','LineWidth',2)
plot(z/(zl/2),(rxz12)/(sqrt(2)*wo),'r','LineWidth',2)
plot(z/(zl/2),(rxz21)/(sqrt(2)*wo),'c','LineWidth',2)
plot(z/(zl/2),(rxz22)/(sqrt(2)*wo),'c','LineWidth',2)
plot(z/(zl/2),-(rxz11)/(sqrt(2)*wo),'r','LineWidth',2)
plot(z/(zl/2),-(rxz12)/(sqrt(2)*wo),'r','LineWidth',2)
plot(z/(zl/2),-(rxz21)/(sqrt(2)*wo),'c','LineWidth',2)
plot(z/(zl/2),-(rxz22)/(sqrt(2)*wo),'c','LineWidth',2)
plot(z/(zl/2), sigmaHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaHx/(sqrt(2)*wo),'b','LineWidth',1.5)
axis1=gca;
set(axis1,'FontSize',28);
xlabel('$z$','Interpreter','latex','FontSize',28) 
ylabel('$x$','Interpreter','latex','FontSize',28)
%graficando la cintura del Elegant Laguerre
hold off
pbaspect([2.5 1 2])

