%Hermite analiticos
clear all
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=10;mu=10;
% Parametros fisicos [micras]
wo=2;
lamb=0.6328;
k=2*pi/lamb;
zl=k*wo^2/2;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl*1.2;%/5;   % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^6+1;        % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=(0:Nz-1)*dz;   % Vector z de propagacion
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
%--------------------Cintura de Hermite en z=0----------------------------%
% Matrices de ceros para las cinturas en x,y en z=0
sigmaHx=zeros(1,sizezp(2));        sigmaHy=zeros(1,sizezp(2));
sigmaH=zeros(1,sizezp(2));
% cinturas en z=0
sigmaHx(1)=sqrt(2*nu+1)*wo/sqrt(2);    sigmaHy(1)=sqrt(2*mu+1)*wo/sqrt(2);
% cintura radial
sigmaH(1)=sqrt(sigmaHx(1)^2+sigmaHy(1));
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
%lineas de la cintura en x,y
line([ sigmaHx(1)/(sqrt(2)*wo), sigmaHx(1)/(sqrt(2)*wo)], [-sigmaHy(1)/(sqrt(2)*wo), sigmaHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',2)
line([-sigmaHx(1)/(sqrt(2)*wo),-sigmaHx(1)/(sqrt(2)*wo)], [-sigmaHy(1)/(sqrt(2)*wo), sigmaHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',2)
line([ sigmaHx(1)/(sqrt(2)*wo),-sigmaHx(1)/(sqrt(2)*wo)], [ sigmaHy(1)/(sqrt(2)*wo), sigmaHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',2)
line([-sigmaHx(1)/(sqrt(2)*wo), sigmaHx(1)/(sqrt(2)*wo)], [-sigmaHy(1)/(sqrt(2)*wo),-sigmaHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',2)
hold off
%-----------Pendiente en x=sigmaeHx(1),y=sigmaeHx(1),z=0------------------%
%-----------------------eHermite para x,y,z-------------------------------%
%-----------Elegant Hermite para x,z=cte y todo y-------------------------%
xi=sigmaHx(1);     yi=x;    zi=0;
% Parametros de Hermite
wz=wo.*sqrt(1+(zi.^2)/(zl^2));
Rz=zi+(zl^2./zi);
Phizx=nu.*atan(zi./zl);
Phizy=mu.*atan(zi./zl);
A=wo./wz;
[HGxzx,NHGxzx]=hermiteg(nu,sqrt(2)*xi./wz);
[HGxzy,NHGxzy]=hermiteg(mu,sqrt(2)*yi./wz);
HGxzx= A.* HGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
NHGxzx=A.*NHGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
HGxzy= A.* HGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
NHGxzy=A.*NHGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
%Funciones Hankel de Hermite 1D para x,z=cte y todo y%
Hxz1x=HGxzx+1i*NHGxzx;
Hxz2x=HGxzx-1i*NHGxzx;
Hxz1y=HGxzy+1i*NHGxzy;
Hxz2y=HGxzy-1i*NHGxzy;
%Funciones Hankel de  Hermite 2D para x,z=cte%
H11xz=Hxz1x.*Hxz1y;
H12xz=Hxz1x.*Hxz2y;
H21xz=Hxz2x.*Hxz1y;
H22xz=Hxz2x.*Hxz2y;
%-----------Elegant Hermite para y,z=cte y todo x-------------------------%
xi=x;     yi=sigmaHy(1);    zi=0;
% Parametros de Hermite
wz=wo.*sqrt(1+(zi.^2)/(zl^2));
Rz=zi+(zl^2./zi);
Phizx=nu.*atan(zi./zl);
Phizy=mu.*atan(zi./zl);
A=wo./wz;
[HGyzx,NHGyzx]=hermiteg(nu,sqrt(2)*xi./wz);
[HGyzy,NHGyzy]=hermiteg(mu,sqrt(2)*yi./wz);
HGyzx= A.* HGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
NHGyzx=A.*NHGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
HGyzy= A.* HGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
NHGyzy=A.*NHGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
%Funciones Hankel de elegant Hermite 1D para y,z=cte%
Hyz1x=HGyzx+1i*NHGyzx;
Hyz2x=HGyzx-1i*NHGyzx;
Hyz1y=HGyzy+1i*NHGyzy;
Hyz2y=HGyzy-1i*NHGyzy;
%Funciones Hankel de elegant Hermite 2D para y,z=cte%
H11yz=Hyz1x.*Hyz1y;
H12yz=Hyz1x.*Hyz2y;
H21yz=Hyz2x.*Hyz1y;
H22yz=Hyz2x.*Hyz2y;
%-----------Elegant Hermite para x,y=cte y todo z-------------------------%
xi=sigmaHx(1); yi=sigmaHy(1);   zi=z;
% Parametros de eHermite
wz=wo.*sqrt(1+(zi.^2)/(zl^2));
Rz=zi+(zl^2./zi);
Phizx=nu.*atan(zi./zl);
Phizy=mu.*atan(zi./zl);
A=wo./wz;
[HGxz,NHGxz]=hermiteg(nu,sqrt(2)*xi./wz);
[HGyz,NHGyz]=hermiteg(mu,sqrt(2)*yi./wz);
HGxz=A.* HGxz .*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
NHGxz=A.*NHGxz.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
HGyz=A.* HGyz .*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
NHGyz=A.*NHGyz.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
%Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
Hx1z=HGxz+1i*NHGxz;
Hx2z=HGxz-1i*NHGxz;
Hy1z=HGyz+1i*NHGyz;
Hy2z=HGyz-1i*NHGyz;
%Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
H11xy=Hx1z.*Hy1z;
H12xy=Hx1z.*Hy2z;
H21xy=Hx2z.*Hy1z;
H22xy=Hx2z.*Hy2z;
% Calculando la pendiente con las derivadas de las fases
% fase x,z constante y todo y
ph11xz=unwrap(angle(H11xz));
ph12xz=unwrap(angle(H12xz));
ph21xz=unwrap(angle(H21xz));
ph22xz=unwrap(angle(H22xz));
% fase y,z constante y todo x
ph11yz=unwrap(angle(H11yz));
ph12yz=unwrap(angle(H12yz));
ph21yz=unwrap(angle(H21yz));
ph22yz=unwrap(angle(H22yz));
% fase x,y constante y todo z
ph11xy=unwrap(angle(H11xy));
ph12xy=unwrap(angle(H12xy));
ph21xy=unwrap(angle(H21xy));
ph22xy=unwrap(angle(H22xy));
% derivada para x,z constante 
wavefront11y=gradient(ph11xz)/dx;
wavefront12y=gradient(ph12xz)/dx;
wavefront21y=gradient(ph21xz)/dx;
wavefront22y=gradient(ph22xz)/dx;
% derivada para y,z constante 
wavefront11x=gradient(ph11yz)/dx;
wavefront12x=gradient(ph12yz)/dx;
wavefront21x=gradient(ph21yz)/dx;
wavefront22x=gradient(ph22yz)/dx;
% derivada para x,y constante
wavefront11z=gradient(ph11xy)/dz+k;
wavefront12z=gradient(ph12xy)/dz+k;
wavefront21z=gradient(ph21xy)/dz+k;
wavefront22z=gradient(ph22xy)/dz+k;

%calculo de las pendietes z con x
m11zxo=wavefront11z(floor(z(1)/dz+1))/wavefront11x(N/2+1+floor(sigmaHx(1)/dx));
m12zxo=wavefront12z(floor(z(1)/dz+1))/wavefront12x(N/2+1+floor(sigmaHx(1)/dx));
m21zxo=wavefront21z(floor(z(1)/dz+1))/wavefront21x(N/2+1+floor(sigmaHx(1)/dx));
m22zxo=wavefront22z(floor(z(1)/dz+1))/wavefront22x(N/2+1+floor(sigmaHx(1)/dx));

%calculo de las pendietes z con y
m11zyo=wavefront11z(floor(z(1)/dz+1))/wavefront11y(N/2+1+floor(sigmaHy(1)/dx)); 
m12zyo=wavefront12z(floor(z(1)/dz+1))/wavefront12y(N/2+1+floor(sigmaHy(1)/dx)); 
m21zyo=wavefront21z(floor(z(1)/dz+1))/wavefront21y(N/2+1+floor(sigmaHy(1)/dx)); 
m22zyo=wavefront22z(floor(z(1)/dz+1))/wavefront22y(N/2+1+floor(sigmaHy(1)/dx)); 

% Vectores para los rayos dos por cada punto en la cintura de cada uno de
% los ejes
% rayos para x
rxz11=zeros(1,sizezp(2));
rxz12=zeros(1,sizezp(2));
rxz21=zeros(1,sizezp(2));
rxz22=zeros(1,sizezp(2));
% rayos para y
ryz11=zeros(1,sizezp(2));
ryz12=zeros(1,sizezp(2));   
ryz21=zeros(1,sizezp(2));  
ryz22=zeros(1,sizezp(2));   

% primer posicion del rayo en x
rx11=sigmaHx(1);
rx12=sigmaHx(1);
rx21=sigmaHx(1);
rx22=sigmaHx(1);

% primer posicion del rayo en y
ry11=sigmaHy(1);
ry12=sigmaHy(1);
ry21=sigmaHy(1);
ry22=sigmaHy(1);

% pendientes en z=0 en el punto(sigmaHx(1),sigmaHy(1))

m11zx=m11zxo;
m12zx=m12zxo;
m21zx=m21zxo;
m22zx=m22zxo;

m11zy=m11zyo;
m12zy=m12zyo;
m21zy=m21zyo;
m22zy=m22zyo;

% rayos para x
rxz11(1)=rx11;
rxz12(1)=rx12;
rxz21(1)=rx21;
rxz22(1)=rx22;
% rayos para y
ryz11(1)=ry11;
ryz12(1)=ry12;   
ryz21(1)=ry21;  
ryz22(1)=ry22;   


vidObj1 = VideoWriter('H.avi');
vidObj1.Quality = 100;
vidObj1.FrameRate = 10;
open(vidObj1);


for ii=2:sizezp(2) %corriendo todos los valores de sp

    % para x
    rx11=(1/m11zx)*(z(ii)-z(ii-1))+rx11;
    rx12=(1/m12zx)*(z(ii)-z(ii-1))+rx12;
    rx21=(1/m21zx)*(z(ii)-z(ii-1))+rx21;
    rx22=(1/m22zx)*(z(ii)-z(ii-1))+rx22;
    % para y
    ry11=(1/m11zy)*(z(ii)-z(ii-1))+ry11;
    ry12=(1/m12zy)*(z(ii)-z(ii-1))+ry12;
    ry21=(1/m21zy)*(z(ii)-z(ii-1))+ry21;
    ry22=(1/m22zy)*(z(ii)-z(ii-1))+ry22;
    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta
    % para x
    rxz11(ii)=rx11;
    rxz12(ii)=rx12;
    rxz21(ii)=rx21;
    rxz22(ii)=rx22;
    % para y
    ryz11(ii)=ry11;
    ryz12(ii)=ry12;
    ryz21(ii)=ry21;
    ryz22(ii)=ry22;
    % Las dos soluciones del Elegant Hermite para x e y para z=cte
    xi=x;
    yi=x;
    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGzx,NHGzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGzy,NHGzy]=hermiteg(mu,sqrt(2)*yi./wz);
    HGzx= A.*HGzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    NHGzx=A.*NHGzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGzy= A.*HGzy .*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);
    NHGzy=A.*NHGzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);
    g=(HGzy')*HGzx;
    gx(:,ii)=g(:,N/2+1);
    gy(:,ii)=g(N/2+1,:);
    %-----------------Pendiente en x=rx11,y=ry11,z=z(ii)------------------%
    %-----------------------eHermite para x,y,z---------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y---------------------%
    xi=rx11;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxzx,NHGxzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGxzy,NHGxzy]=hermiteg(mu,sqrt(2)*yi./wz);
    HGxzx= A.*HGxzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGxzx=A.*NHGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGxzy= A.*HGxzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
    NHGxzy=A.*NHGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz1x=HGxzx+1i*NHGxzx;
    Hxz1y=HGxzy+1i*NHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H11xz=Hxz1x.*Hxz1y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry11;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGyzx,NHGyzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyzy,NHGyzy]=hermiteg(mu,sqrt(2)*yi./wz);
    HGyzx= A.*HGyzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    NHGyzx=A.*NHGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGyzy= A.*HGyzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    NHGyzy=A.*NHGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz1x=HGyzx+1i*NHGyzx;
    Hyz1y=HGyzy+1i*NHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H11yz=Hyz1x.*Hyz1y;
    %-----------Elegant Hermite para x,y=cte y todo z----------------------%
    xi=rx11; yi=ry11;   zi=z;
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxz,NHGxz]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyz,NHGyz]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxz= A.*HGxz.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGxz=A.*NHGxz.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGyz= A.*HGyz.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy); 
    NHGyz=A.*NHGyz.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx1z=HGxz+1i*NHGxz;
    Hy1z=HGyz+1i*NHGyz;
    %Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
    H11xy=Hx1z.*Hy1z;
    % Calculando la pendiente con las derivadas de las fases
    % fase x,z constante y todo y
    ph11xz=unwrap(angle(H11xz));
    % fase y,z constante y todo x
    ph11yz=unwrap(angle(H11yz));
    % fase x,y constante y todo z
    ph11xy=unwrap(angle(H11xy));
    % derivada para x,z constante 
    wavefront11y=gradient(ph11xz)/dx;
    % derivada para y,z constante 
    wavefront11x=gradient(ph11yz)/dx;
    % derivada para x,y constante
    wavefront11z=gradient(ph11xy)/dz+k;

    %calculo de la pendiete z con x
    m11zx=wavefront11z(floor(z(ii)/dz+1))/wavefront11x(N/2+1+floor(rx11/dx));

    %calculo de la pendiete z con y
    m11zy=wavefront11z(floor(z(ii)/dz+1))/wavefront11y(N/2+1+floor(ry11/dx)); 
    %------------------Pendiente en x=rx12,y=ry12,z=z(ii)-----------------%
    %-----------------------eHermite para x,y,z---------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y---------------------%
    xi=rx12;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxzx,NHGxzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGxzy,NHGxzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxzx= A.*HGxzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGxzx=A.*NHGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGxzy= A.*HGxzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    NHGxzy=A.*NHGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz1x=HGxzx+1i*NHGxzx;
    Hxz2y=HGxzy-1i*NHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H12xz=Hxz1x.*Hxz2y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry12;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGyzx,NHGyzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyzy,NHGyzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGyzx=A .*HGyzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGyzx=A.*NHGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGyzy=A .*HGyzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyzy=A.*NHGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz1x=HGyzx+1i*NHGyzx;
    Hyz2y=HGyzy-1i*NHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H12yz=Hyz1x.*Hyz2y;
    %-----------Elegant Hermite para x,y=cte y todo z-------------------------%
    xi=rx12; yi=ry12;   zi=z;
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxz,NHGxz]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyz,NHGyz]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxz= A.*HGxz.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);  
    NHGxz=A.*NHGxz.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    HGyz= A.*HGyz.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyz=A.*NHGyz.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx1z=HGxz+i*NHGxz;
    Hy2z=HGyz-i*NHGyz;
    %Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
    H12xy=Hx1z.*Hy2z;
    % Calculando la pendiente con las derivadas de las fases
    % fase x,z constante y todo y
    ph12xz=unwrap(angle(H12xz));
    % fase y,z constante y todo x
    ph12yz=unwrap(angle(H12yz));
    % fase x,y constante y todo z
    ph12xy=unwrap(angle(H12xy));
    % derivada para x,z constante 
    wavefront12y=gradient(ph12xz)/dx;
    % derivada para y,z constante 
    wavefront12x=gradient(ph12yz)/dx;
    % derivada para x,y constante
    wavefront12z=gradient(ph12xy)/dz+k;

    %calculo de las pendietes z con x
    m12zx=wavefront12z(floor(z(ii)/dz+1))/wavefront12x(N/2+1+floor(rx12/dx));

    %calculo de las pendietes z con y
    m12zy=wavefront12z(floor(z(ii)/dz+1))/wavefront12y(N/2+1+floor(ry12/dx)); 

    %-------------------Pendiente en x=rx21,y=rx21,z=z(ii)--------------------%
    %-----------------------eHermite para x,y,z---------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y---------------------%
    xi=rx21;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxzx,NHGxzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGxzy,NHGxzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxzx= A.*HGxzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGxzx=A.*NHGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGxzy= A.*HGxzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    NHGxzy=A.*NHGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz2x=HGxzx-1i*NHGxzx;
    Hxz1y=HGxzy+1i*NHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H21xz=Hxz2x.*Hxz1y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry21;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGyzx,NHGyzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyzy,NHGyzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGyzx=A .*HGyzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGyzx=A.*NHGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGyzy=A .*HGyzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyzy=A.*NHGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz2x=HGyzx-i*NHGyzx;
    Hyz1y=HGyzy+i*NHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H21yz=Hyz2x.*Hyz1y;
    %-----------Elegant Hermite para x,y=cte y todo z-------------------------%
    xi=rx21; yi=ry21;   zi=z;
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxz,NHGxz]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyz,NHGyz]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxz= A.*HGxz.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);  
    NHGxz=A.*NHGxz.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    HGyz= A.*HGyz.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyz=A.*NHGyz.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx2z=HGxz-i*NHGxz;
    Hy1z=HGyz+i*NHGyz;
    %Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
    H21xy=Hx2z.*Hy1z;
    % Calculando la pendiente con las derivadas de las fases
    % fase x,z constante y todo y
    ph21xz=unwrap(angle(H21xz));
    % fase y,z constante y todo x
    ph21yz=unwrap(angle(H21yz));
    % fase x,y constante y todo z
    ph21xy=unwrap(angle(H21xy));
    % derivada para x,z constante 
    wavefront21y=gradient(ph21xz)/dx;
    % derivada para y,z constante 
    wavefront21x=gradient(ph21yz)/dx;
    % derivada para x,y constante
    wavefront21z=gradient(ph21xy)/dz+k;

    %calculo de las pendietes z con x
    m21zx=wavefront21z(floor(z(ii)/dz+1))/wavefront21x(N/2+1+floor(rx21/dx));

    %calculo de las pendietes z con y
    m21zy=wavefront21z(floor(z(ii)/dz+1))/wavefront21y(N/2+1+floor(ry21/dx)); 

    %-------------------Pendiente en x=rx22,y=ry22,z=z(ii)----------------%
    %-----------------------eHermite para x,y,z---------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y---------------------%
    xi=rx22;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxzx,NHGxzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGxzy,NHGxzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxzx= A.*HGxzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGxzx=A.*NHGxzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGxzy= A.*HGxzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    NHGxzy=A.*NHGxzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz2x=HGxzx-1i*NHGxzx;
    Hxz2y=HGxzy-1i*NHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H22xz=Hxz2x.*Hxz2y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry22;    zi=z(ii);
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGyzx,NHGyzx]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyzy,NHGyzy]=hermiteg(nu,sqrt(2)*yi./wz);
    HGyzx=A .*HGyzx.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    NHGyzx=A.*NHGyzx.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);
    HGyzy=A .*HGyzy.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyzy=A.*NHGyzy.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz2x=HGyzx-i*NHGyzx;
    Hyz2y=HGyzy-i*NHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H22yz=Hyz2x.*Hyz2y;
    %-----------Elegant Hermite para x,y=cte y todo z---------------------%
    xi=rx22; yi=ry22;   zi=z;
    % Parametros de eHermite
    wz=wo.*sqrt(1+(zi.^2)/(zl^2));
    Rz=zi+(zl^2./zi);
    Phizx=nu.*atan(zi./zl);
    Phizy=mu.*atan(zi./zl);
    A=wo./wz;
    [HGxz,NHGxz]=hermiteg(nu,sqrt(2)*xi./wz);
    [HGyz,NHGyz]=hermiteg(nu,sqrt(2)*yi./wz);
    HGxz= A.*HGxz.* exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx);  
    NHGxz=A.*NHGxz.*exp(-1i*k*xi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizx); 
    HGyz= A.*HGyz.* exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);   
    NHGyz=A.*NHGyz.*exp(-1i*k*yi.^2./(2*Rz)).*exp(-1i*k*zi).*exp(1i.*k.*Phizy);  
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx2z=HGxz-i*NHGxz;
    Hy2z=HGyz-i*NHGyz;
    %Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
    H22xy=Hx2z.*Hy2z;
    % Calculando la pendiente con las derivadas de las fases
    % fase x,z constante y todo y
    ph22xz=unwrap(angle(H22xz));
    % fase y,z constante y todo x
    ph22yz=unwrap(angle(H22yz));
    % fase x,y constante y todo z
    ph22xy=unwrap(angle(H22xy));
    % derivada para x,z constante 
    wavefront22y=gradient(ph22xz)/dx;
    % derivada para y,z constante 
    wavefront22x=gradient(ph22yz)/dx;
    % derivada para x,y constante
    wavefront22z=gradient(ph22xy)/dz+k;

    %calculo de las pendietes z con x
    m22zx=wavefront22z(floor(z(ii)/dz+1))/wavefront22x(N/2+1+floor(rx22/dx));

    %calculo de las pendietes z con y
    m22zy=wavefront22z(floor(z(ii)/dz+1))/wavefront22y(N/2+1+floor(ry22/dx)); 

%     figure(1)
%     plot(ph11xz)
%     figure(2)
%     plot(ph12xz)
%     figure(3)
%     plot(ph21xz)
%     figure(4)
%     plot(ph22xz)
    figure(10)
    fig10=figure(10);
    fig10.Name=([' z = ',num2str(z(ii)/(sqrt(2)*wo))]);
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs((HGzy')*(HGzx)))
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis1=gca;
    set(axis1,'FontSize',28);
    xlabel('$x$','Interpreter','latex','FontSize',28) 
    ylabel('$y$','Interpreter','latex','FontSize',28)
    hold on
    %lineas de la cintura en x,y
    line([ sigmaHx(ii)/(sqrt(2)*wo), sigmaHx(ii)/(sqrt(2)*wo)], [-sigmaHy(ii)/(sqrt(2)*wo), sigmaHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([-sigmaHx(ii)/(sqrt(2)*wo),-sigmaHx(ii)/(sqrt(2)*wo)], [-sigmaHy(ii)/(sqrt(2)*wo), sigmaHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([ sigmaHx(ii)/(sqrt(2)*wo),-sigmaHx(ii)/(sqrt(2)*wo)], [ sigmaHy(ii)/(sqrt(2)*wo), sigmaHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([-sigmaHx(ii)/(sqrt(2)*wo), sigmaHx(ii)/(sqrt(2)*wo)], [-sigmaHy(ii)/(sqrt(2)*wo),-sigmaHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
%     line([ sigmaeHx(1)/(sqrt(2)*wo), sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
%     line([-sigmaeHx(1)/(sqrt(2)*wo),-sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
%     line([ sigmaeHx(1)/(sqrt(2)*wo),-sigmaeHx(1)/(sqrt(2)*wo)], [ sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
%     line([-sigmaeHx(1)/(sqrt(2)*wo), sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo),-sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
%     % Dibujando circulo de valor esperado
%     re=vesperado(ii)/(sqrt(2)*wo);
%     xp=re*cos(ang);
%     yp=re*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','m');
%     % Valores esperados en x, y
%     line([ vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([-vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([ vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [ vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([-vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo),-vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    %esquina superior derecha
    plot(rx11/(sqrt(2)*wo),ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(rx12/(sqrt(2)*wo),ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx21/(sqrt(2)*wo),ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx22/(sqrt(2)*wo),ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %esquina superior izquierda
    plot(-rx11/(sqrt(2)*wo),ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(-rx12/(sqrt(2)*wo),ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx21/(sqrt(2)*wo),ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx22/(sqrt(2)*wo),ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %esquina inferior izquierna
    plot(-rx11/(sqrt(2)*wo),-ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(-rx12/(sqrt(2)*wo),-ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx21/(sqrt(2)*wo),-ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(-rx22/(sqrt(2)*wo),-ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    %esquina inferior derecha
    plot(rx11/(sqrt(2)*wo),-ry11/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','r')
    plot(rx12/(sqrt(2)*wo),-ry12/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx21/(sqrt(2)*wo),-ry21/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','m')
    plot(rx22/(sqrt(2)*wo),-ry22/(sqrt(2)*wo),'+','MarkerSize',15,'LineWidth',3,'color','c')
    hold off
    pause(.001)
    writeVideo(vidObj1, getframe(gca));
end  
close(vidObj1);
figure(8)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gx).^.75)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en y')
hold on
%graficando las rectas de las normales al frente de onda
%graficando los rayos de las normales al frente de onda
plot(z/(zl/2),(rxz11)/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),(rxz12)/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),(rxz21)/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),(rxz22)/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),-(rxz11)/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),-(rxz12)/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),-(rxz21)/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),-(rxz22)/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2), sigmaHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaHx/(sqrt(2)*wo),'b','LineWidth',1.5)
    axis1=gca;
    set(axis1,'FontSize',28);
    xlabel('$z$','Interpreter','latex','FontSize',28) 
    ylabel('$x$','Interpreter','latex','FontSize',28)
%graficando la cintura del Elegant Laguerre
hold off
pbaspect([2.5 1 2])

figure(9)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.15)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en y')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),((1/m11zyo)*z+sigmaHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m12zyo)*z+sigmaHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m21zyo)*z+sigmaHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),((1/m22zyo)*z+sigmaHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%graficando los rayos de las normales al frente de onda
plot(z/(zl/2),(ryz11)/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),(ryz12)/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),(ryz21)/(sqrt(2)*wo),'m','LineWidth',1.5)
plot(z/(zl/2),(ryz22)/(sqrt(2)*wo),'m','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
hold off
pbaspect([2.5 1 2])

