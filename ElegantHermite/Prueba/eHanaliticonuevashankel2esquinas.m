%eHermite analiticos
clear all
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=10;mu=20;
% Parametros fisicos [micras]
wo=2;
lamb=0.6328;    %micras
k=2*pi/lamb;
zl=k*wo^2/2;
qo=-1i*zl;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % N�mero de puntos para x, y
n=-N/2:N/2-1;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl*1.2;%/5; % Tama�o de la ventana en z (distancia a la cual propagar)
Nz=2^6+1;          % N�mero de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=(0:Nz-1)*dz;       % Vector z de propagacion
% calculando el tama�o de zp
sizezp=size(z);
% Cintura de gaussiana a la distancia del tama�o de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2);
% Cintura de Laguerre a la distancia del tama�o de la ventana en z
niu=max(nu,mu);
sigmaH=ws*sqrt((2*niu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tama�o de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(sigmaH)*2;%*1.6/(1.97)^2;  % Tama�o de la ventana
dx=Dx/N;            % Resoluci�n
x=n*dx;             % Vector x (plano de campo a propapar)
y=x;
[X]=meshgrid(x);
x1=x(N/2:N);    % vector de 0 a Dx
% Muestreo del vector de frecuencia
Du=1/dx;        % Tama�o de la ventana de Fourier
du=1/Dx;        % Resoluci�n
u=n*du;         % Vector
[U]=meshgrid(u);
% Vectores kx,ky
kx=2*pi*u;
[Kx]=meshgrid(kx); 
% Para dibujar circulos 
ang=0:0.01:2*pi;
%----------------Cintura de elegant Hermite en z=0------------------------%
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=x;    yi=x;    zi=0;
% Parametros de eHermite
An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
[eHGzx,eNHGzx]=ehermite(nu,(alpha.*xi));
[eHGzy,eNHGzy]=ehermite(mu,(alpha.*yi));
eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy; 
g=(eHGzy')*eHGzx;
% Matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));
% Guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);    gy(:,1)=g(:,N/2+1);

% Matrices de ceros para las cinturas en x,y en z=0
sigmaeHx=zeros(1,sizezp(2));        sigmaeHy=zeros(1,sizezp(2));
sigmaeH=zeros(1,sizezp(2));
% Intensidades en total y en x, y en z=0
% Ix=eHG(N/2+1,:).*conj(eHG(N/2+1,:)); Iy=eHG(:,N/2+1).*conj(eHG(:,N/2+1));
Ix=eHGzx.*conj(eHGzx);    Iy=eHGzy.*conj(eHGzy);
I=(Iy)'*Ix;
% Calculando las cinturas
sigmaeHx(1)=sqrt(2*(trapz(x,Ix.*(x.^2)))./(trapz(x,Ix)));
sigmaeHy(1)=sqrt(2*(trapz(x,Iy.*(x.^2)))./(trapz(x,Iy)));
% cintura radial
sigmaeH(1)=sqrt(2*(trapz(y,trapz(x,I.*(X.^2+(X').^2),2)))./(trapz(y,trapz(x,I,2))));
%--------------------------campo inicial----------------------------------%
figure(1)
imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^1.1)
axis ([-2.5,2.5,-2.5,2.5])
axis square
colormap(mapgreen)
set(gca,'YDir','normal')
axis1=gca;
set(axis1,'FontSize',28);
xlabel('$x$','Interpreter','latex','FontSize',28) 
ylabel('$y$','Interpreter','latex','FontSize',28)
hold on
%lineas de la cintura en x,y
line([ sigmaeHx(1)/(sqrt(2)*wo), sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
line([-sigmaeHx(1)/(sqrt(2)*wo),-sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
line([ sigmaeHx(1)/(sqrt(2)*wo),-sigmaeHx(1)/(sqrt(2)*wo)], [ sigmaeHy(1)/(sqrt(2)*wo), sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
line([-sigmaeHx(1)/(sqrt(2)*wo), sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo),-sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
hold off


H1x=eHGzx+i*eNHGzx;
H1y=eHGzx+i*eNHGzx;

H1=(H1y')*(H1x);
%fase
phH = phase_unwrap(angle(H1));
figure(2)
mesh(x/(sqrt(2)*wo),x/(sqrt(2)*wo),phH)
set(gca,'YDir','normal')
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',25) 
ylabel('$y$','Interpreter','latex','FontSize',25)
zlabel('$phase$','Interpreter','latex','FontSize',25)
axis([-5 5 -5 5 -35 35])
%------------------Calculo del valor esperado en z=0----------------------%
% Matriz de ceros para el valor esperado
vesperado=zeros(1,sizezp(2));
vesperadox=zeros(1,sizezp(2));
vesperadoy=zeros(1,sizezp(2));
% Calculando el valor esperado
vesperado(1)=((trapz(y,trapz(x,I.*sqrt(X.^2+(X').^2),2)))./(trapz(y,trapz(x,I,2))));
% Tomando intensidad de 0 a Dx y de 0 a Dy (solo un cuadrante)  
Ix1=Ix(N/2:N);
Iy1=Iy(N/2:N);
% Valor esperado en x,y de un cuadrante
vesperadox(1)=(trapz(x1,Ix1.*x1))./(trapz(x1,Ix1));
vesperadoy(1)=(trapz(x1,(Iy1).*x1))./(trapz(x1,Iy1));
%-------------------Pendiente en x=sigmaeHx(1),y=0,z=0--------------------%
%-----------------------eHermite para x,y,z-------------------------------%
%-----------Elegant Hermite para x,z=cte y todo y-------------------------%
xi=sigmaeHx(1);     yi=x;    zi=0;
% Parametros de eHermite
qz=zi+1i*zl;    C=1;
An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
[eHGxzx,eNHGxzx]=ehermite(nu,(alpha.*xi));
[eHGxzy,eNHGxzy]=ehermite(mu,(alpha.*yi));
eHGxzx=An.*eHGxzx; eNHGxzx=An.*eNHGxzx;
eHGxzy=Am.*eHGxzy; eNHGxzy=Am.*eNHGxzy; 
%Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
Hxz1x=eHGxzx+1i*eNHGxzx;
Hxz2x=eHGxzx-1i*eNHGxzx;
Hxz1y=eHGxzy+1i*eNHGxzy;
Hxz2y=eHGxzy-1i*eNHGxzy;
%Funciones Hankel de elegant Hermite 2D para x,z=cte%
H11xz=Hxz1x.*Hxz1y;
H12xz=Hxz1x.*Hxz2y;
H21xz=Hxz2x.*Hxz1y;
H22xz=Hxz2x.*Hxz2y;
%-----------Elegant Hermite para y,z=cte y todo x-------------------------%
xi=x;     yi=0;    zi=0;
% Parametros de eHermite
qz=zi+i*zl;    C=1;
An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
[eHGyzx,eNHGyzx]=ehermite(nu,(alpha.*xi));
[eHGyzy,eNHGyzy]=ehermite(mu,(alpha.*yi));
eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGyzx;
eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGyzy; 
%Funciones Hankel de elegant Hermite 1D para y,z=cte%
Hyz1x=eHGyzx+i*eNHGyzx;
Hyz2x=eHGyzx-i*eNHGyzx;
Hyz1y=eHGyzy+i*eNHGyzy;
Hyz2y=eHGyzy-i*eNHGyzy;
%Funciones Hankel de elegant Hermite 2D para y,z=cte%
H11yz=Hyz1x.*Hyz1y;
H12yz=Hyz1x.*Hyz2y;
H21yz=Hyz2x.*Hyz1y;
H22yz=Hyz2x.*Hyz2y;
%-----------Elegant Hermite para x,y=cte y todo z-------------------------%
xi=sigmaeHx(1); yi=0;   zi=z;
% Parametros de eHermite
qz=zi+i*zl; 
An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
[eHGxz,eNHGxz]=ehermite(nu,(alpha.*xi));
[eHGyz,eNHGyz]=ehermite(mu,(alpha.*yi));
eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
%Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
Hx1z=eHGxz+i*eNHGxz;
Hx2z=eHGxz-i*eNHGxz;
Hy1z=eHGyz+i*eNHGyz;
Hy2z=eHGyz-i*eNHGyz;
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
m11zxo=wavefront11z(floor(z(1)/dz+1))/wavefront11x(N/2+1+floor(sigmaeHx(1)/dx));
m12zxo=wavefront12z(floor(z(1)/dz+1))/wavefront12x(N/2+1+floor(sigmaeHx(1)/dx));
m21zxo=wavefront21z(floor(z(1)/dz+1))/wavefront21x(N/2+1+floor(sigmaeHx(1)/dx));
m22zxo=wavefront22z(floor(z(1)/dz+1))/wavefront22x(N/2+1+floor(sigmaeHx(1)/dx));

%calculo de las pendietes z con y
m11zyo=wavefront11z(floor(z(1)/dz+1))/wavefront11y(N/2+1+floor(sigmaeHy(1)/dx)); 
m12zyo=wavefront12z(floor(z(1)/dz+1))/wavefront12y(N/2+1+floor(sigmaeHy(1)/dx)); 
m21zyo=wavefront21z(floor(z(1)/dz+1))/wavefront21y(N/2+1+floor(sigmaeHy(1)/dx)); 
m22zyo=wavefront22z(floor(z(1)/dz+1))/wavefront22y(N/2+1+floor(sigmaeHy(1)/dx)); 


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
rx11=sigmaeHx(1);
rx12=sigmaeHx(1);
rx21=sigmaeHx(1);
rx22=sigmaeHx(1);

% primer posicion del rayo en y
ry11=sigmaeHy(1);
ry12=sigmaeHy(1);
ry21=sigmaeHy(1);
ry22=sigmaeHy(1);

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


vidObj1 = VideoWriter('eH.avi');
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
    % Guardando la posici�n de este rayo para z(ii) para generar la curva
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
    qz=zi+i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGzx,eNHGzx]=ehermite(nu,(alpha.*xi));
    [eHGzy,eNHGzy]=ehermite(mu,(alpha.*yi));
    eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
    eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy;
    g=(eHGzy')*eHGzx;
    gx(:,ii)=g(:,N/2+1);
    gy(:,ii)=g(N/2+1,:);
    %---------------Calculando las cinturas en x,y en z(ii)---------------%
    Iy= eHGzy.*conj(eHGzy);	%intensidad
    Ix= eHGzx.*conj(eHGzx);	%intensidad
    I=Iy'*Ix;
    sigmaeHx(ii)=sqrt(2*(trapz(x,Ix.*(x.^2)))./(trapz(x,Ix)));
    sigmaeHy(ii)=sqrt(2*(trapz(x,Iy.*(x.^2)))./(trapz(x,Iy)));
    %------------------Calculando el valor esperado-----------------------%
    %valor esperado radial
    vesperado(ii)=((trapz(y,trapz(x,I.*sqrt(X.^2+(X').^2),2)))./(trapz(y,trapz(x,I,2))));
    % Intensidad en un cuadrante par calcular el valor esperado en las
    % componentes
    Ix1=Ix(N/2:N);  Iy1=Iy(N/2:N);  %I1=Iy1*Ix1;
    vesperadoy(ii)=(trapz(x1,Ix1.*x1))./(trapz(x1,Ix1));
    vesperadox(ii)=(trapz(x1,(Iy1).*x1))./(trapz(x1,Iy1));
    %-----------------Pendiente en x=rx11,y=ry11,z=z(ii)------------------%
    %-----------------------eHermite para x,y,z---------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y---------------------%
    xi=rx11;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+1i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxzx,eNHGxzx]=ehermite(nu,(alpha.*xi));
    [eHGxzy,eNHGxzy]=ehermite(mu,(alpha.*yi));
    eHGxzx=An.*eHGxzx; eNHGxzx=An.*eNHGxzx;
    eHGxzy=Am.*eHGxzy; eNHGxzy=Am.*eNHGxzy; 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz1x=eHGxzx+1i*eNHGxzx;
    Hxz1y=eHGxzy+1i*eNHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H11xz=Hxz1x.*Hxz1y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry11;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGyzx,eNHGyzx]=ehermite(nu,(alpha.*xi));
    [eHGyzy,eNHGyzy]=ehermite(mu,(alpha.*yi));
    eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGyzx;
    eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGyzy; 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz1x=eHGyzx+i*eNHGyzx;
    Hyz1y=eHGyzy+i*eNHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H11yz=Hyz1x.*Hyz1y;
    %-----------Elegant Hermite para x,y=cte y todo z----------------------%
    xi=rx11; yi=ry11;   zi=z;
    % Parametros de eHermite
    qz=zi+i*zl; 
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxz,eNHGxz]=ehermite(nu,(alpha.*xi));
    [eHGyz,eNHGyz]=ehermite(mu,(alpha.*yi));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx1z=eHGxz+i*eNHGxz;
    Hy1z=eHGyz+i*eNHGyz;
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
    qz=zi+1i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxzx,eNHGxzx]=ehermite(nu,(alpha.*xi));
    [eHGxzy,eNHGxzy]=ehermite(mu,(alpha.*yi));
    eHGxzx=An.*eHGxzx; eNHGxzx=An.*eNHGxzx;
    eHGxzy=Am.*eHGxzy; eNHGxzy=Am.*eNHGxzy; 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz1x=eHGxzx+1i*eNHGxzx;
    Hxz2y=eHGxzy-1i*eNHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H12xz=Hxz1x.*Hxz2y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry12;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGyzx,eNHGyzx]=ehermite(nu,(alpha.*xi));
    [eHGyzy,eNHGyzy]=ehermite(mu,(alpha.*yi));
    eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGyzx;
    eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGyzy; 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz1x=eHGyzx+i*eNHGyzx;
    Hyz2y=eHGyzy-i*eNHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H12yz=Hyz1x.*Hyz2y;
    %-----------Elegant Hermite para x,y=cte y todo z-------------------------%
    xi=rx12; yi=ry12;   zi=z;
    % Parametros de eHermite
    qz=zi+i*zl; 
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxz,eNHGxz]=ehermite(nu,(alpha.*xi));
    [eHGyz,eNHGyz]=ehermite(mu,(alpha.*yi));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx1z=eHGxz+i*eNHGxz;
    Hy2z=eHGyz-i*eNHGyz;
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
    qz=zi+1i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxzx,eNHGxzx]=ehermite(nu,(alpha.*xi));
    [eHGxzy,eNHGxzy]=ehermite(mu,(alpha.*yi));
    eHGxzx=An.*eHGxzx; eNHGxzx=An.*eNHGxzx;
    eHGxzy=Am.*eHGxzy; eNHGxzy=Am.*eNHGxzy; 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz2x=eHGxzx-1i*eNHGxzx;
    Hxz1y=eHGxzy+1i*eNHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H21xz=Hxz2x.*Hxz1y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry21;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGyzx,eNHGyzx]=ehermite(nu,(alpha.*xi));
    [eHGyzy,eNHGyzy]=ehermite(mu,(alpha.*yi));
    eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGyzx;
    eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGyzy; 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz2x=eHGyzx-i*eNHGyzx;
    Hyz1y=eHGyzy+i*eNHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H21yz=Hyz2x.*Hyz1y;
    %-----------Elegant Hermite para x,y=cte y todo z-------------------------%
    xi=rx21; yi=ry21;   zi=z;
    % Parametros de eHermite
    qz=zi+i*zl; 
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxz,eNHGxz]=ehermite(nu,(alpha.*xi));
    [eHGyz,eNHGyz]=ehermite(mu,(alpha.*yi));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx2z=eHGxz-i*eNHGxz;
    Hy1z=eHGyz+i*eNHGyz;
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
    qz=zi+1i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxzx,eNHGxzx]=ehermite(nu,(alpha.*xi));
    [eHGxzy,eNHGxzy]=ehermite(mu,(alpha.*yi));
    eHGxzx=An.*eHGxzx; eNHGxzx=An.*eNHGxzx;
    eHGxzy=Am.*eHGxzy; eNHGxzy=Am.*eNHGxzy; 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz2x=eHGxzx-1i*eNHGxzx;
    Hxz2y=eHGxzy-1i*eNHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    H22xz=Hxz2x.*Hxz2y;
    %-----------Elegant Hermite para y,z=cte y todo x---------------------%
    xi=x;     yi=ry22;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGyzx,eNHGyzx]=ehermite(nu,(alpha.*xi));
    [eHGyzy,eNHGyzy]=ehermite(mu,(alpha.*yi));
    eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGyzx;
    eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGyzy; 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz2x=eHGyzx-i*eNHGyzx;
    Hyz2y=eHGyzy-i*eNHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    H22yz=Hyz2x.*Hyz2y;
    %-----------Elegant Hermite para x,y=cte y todo z---------------------%
    xi=rx22; yi=ry22;   zi=z;
    % Parametros de eHermite
    qz=zi+i*zl; 
    An=(1+zi.^2/zl^2).^(-(nu+1)/4).*(exp(-(i/2)*(nu+1).*atan(zi./zl)));
    Am=(1+zi.^2/zl^2).^(-(mu+1)/4).*(exp(-(i/2)*(mu+1).*atan(zi./zl)));
    alpha=sqrt(k/2)*((zi.^2+zl^2).^(-1/4)).*exp(-(i/2).*atan(zi./zl));
    [eHGxz,eNHGxz]=ehermite(nu,(alpha.*xi));
    [eHGyz,eNHGyz]=ehermite(mu,(alpha.*yi));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx2z=eHGxz-i*eNHGxz;
    Hy2z=eHGyz-i*eNHGyz;
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
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs((eHGzy')*(eHGzx)))
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis1=gca;
    set(axis1,'FontSize',28);
    xlabel('$x$','Interpreter','latex','FontSize',28) 
    ylabel('$y$','Interpreter','latex','FontSize',28)
    hold on
    %lineas de la cintura en x,y
    line([ sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([-sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([ sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [ sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
    line([-sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo),-sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',3)
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
plot(z/(zl/2), sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
    axis1=gca;
    set(axis1,'FontSize',28);
    xlabel('$z$','Interpreter','latex','FontSize',28) 
    ylabel('$x$','Interpreter','latex','FontSize',28)
%graficando la cintura del Elegant Laguerre
hold off
pbaspect([2.5 1 2])

% figure(9)
% pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.15)
% shading interp
% colormap(mapgreen)
% % title('Campo propagado corte en y')
% hold on
% %graficando las rectas de las normales al frente de onda
% plot(z/(zl/2),((1/m11zyo)*z+sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
% plot(z/(zl/2),((1/m12zyo)*z+sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
% plot(z/(zl/2),((1/m21zyo)*z+sigmaeHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
% plot(z/(zl/2),((1/m22zyo)*z+sigmaeHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
% %graficando los rayos de las normales al frente de onda
% plot(z/(zl/2),(ryz11)/(sqrt(2)*wo),'b','LineWidth',1.5)
% plot(z/(zl/2),(ryz12)/(sqrt(2)*wo),'b','LineWidth',1.5)
% plot(z/(zl/2),(ryz21)/(sqrt(2)*wo),'m','LineWidth',1.5)
% plot(z/(zl/2),(ryz22)/(sqrt(2)*wo),'m','LineWidth',1.5)
% %graficando la cintura del Elegant Laguerre
% hold off
% pbaspect([2.5 1 2])
% 
