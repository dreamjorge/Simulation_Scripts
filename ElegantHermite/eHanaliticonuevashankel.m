%eHermite analiticos
clear all
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=20;mu=20;
% Parametros fisicos [micras]
wo=2;
lamb=0.6328;    %micras
k=2*pi/lamb;
zl=k*wo^2/2;
qo=1i*zl;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl;%/5; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^8;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=(0:Nz-1)*dz;       % Vector z de propagacion
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
Dx=(sigmaH)*2.7;%*1.6/(1.97)^2;  % Tamaño de la ventana
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
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=x;    yi=x;    zi=0;
% Parametros de eHermite
qz=zi+1i*zl;    C=1;
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
alphax=sqrt(1i*(k)./(2*qz)).*xi;
alphay=sqrt(1i*(k)./(2*qz)).*yi;
[eHGzx,eNHGzx]=ehermite(nu,(alphax));
[eHGzy,eNHGzy]=ehermite(mu,(alphay));
eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy; 
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
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
alphax=sqrt(1i*(k)./(2*qz)).*xi;
alphay=sqrt(1i*(k)./(2*qz)).*yi;
[eHGxzx,eNHGxzx]=ehermite(nu,(alphax));
[eHGxzy,eNHGxzy]=ehermite(mu,(alphay));
eHGxzx=An.*eHGzx; eNHGxzx=An.*eNHGxzx;
eHGxzy=Am.*eHGzy; eNHGxzy=Am.*eNHGxzy; 
%Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
Hxz1x=eHGxzx+1i*eNHGxzx;
Hxz2x=eHGxzx-1i*eNHGxzx;
Hxz1y=eHGxzy+1i*eNHGxzy;
Hxz2y=eHGxzy-1i*eNHGxzy;
%Funciones Hankel de elegant Hermite 2D para x,z=cte%
Hxz1yx=(eHGxzy).*Hxz1x;
Hxz2yx=(eHGxzy).*Hxz2x;
Hxz1xy= (Hxz1y).*eHGxzx;
Hxz2xy= (Hxz2y).*eHGxzx;
%-----------Elegant Hermite para y,z=cte y todo x-------------------------%
xi=x;     yi=sigmaeHy(1);    zi=0;
% Parametros de eHermite
qz=zi+i*zl;    C=1;
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
alphax=sqrt(i*(k)./(2*qz)).*xi;
alphay=sqrt(i*(k)./(2*qz)).*yi;
[eHGyzx,eNHGyzx]=ehermite(nu,(alphax));
[eHGyzy,eNHGyzy]=ehermite(mu,(alphay));
eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGzx;
eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGzy; 
%Funciones Hankel de elegant Hermite 1D para y,z=cte%
Hyz1x=eHGyzx+i*eNHGyzx;
Hyz2x=eHGyzx-i*eNHGyzx;
Hyz1y=eHGyzy+i*eNHGyzy;
Hyz2y=eHGyzy-i*eNHGyzy;
%Funciones Hankel de elegant Hermite 2D para y,z=cte%
Hyz1yx=(eHGyzy).*Hyz1x;
Hyz2yx=(eHGyzy).*Hyz2x;
Hyz1xy=(Hyz1y).*eHGyzx;
Hyz2xy=(Hyz2y).*eHGyzx;
%-----------Elegant Hermite para x,y=cte y todo z-------------------------%
xi=sigmaeHx(1); yi=sigmaeHy(1);   zi=z;
% Parametros de eHermite
qz=zi+i*zl; 
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
alphax=sqrt(i*(k)./(2*qz)).*xi;
alphay=sqrt(i*(k)./(2*qz)).*yi;
[eHGxz,eNHGxz]=ehermite(nu,(alphax));
[eHGyz,eNHGyz]=ehermite(mu,(alphay));
eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
%Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
Hx1z=eHGxz+i*eNHGxz;
Hx2z=eHGxz-i*eNHGxz;
Hy1z=eHGyz+i*eNHGyz;
Hy2z=eHGyz-i*eNHGyz;
%Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
Hyx1z=(eHGyz).*Hx1z;
Hyx2z=(eHGyz).*Hx2z;
Hxy1z=(Hy1z).*eHGxz;
Hxy2z=(Hy2z).*eHGxz;
% Calculando la pendiente con las derivadas de las fases
% fase x,z constante y todo y
phxzy=unwrap(angle(Hxz2xy));
% fase y,z constante y todo x
phyzx=unwrap(angle(Hyz2yx));
% fase x,y constante y todo z
phxyz=unwrap(angle(Hxy2z));
% derivada para x,z constante 
wavefronty=gradient(phxzy)/dx;
% derivada para y,z constante 
wavefrontx=gradient(phyzx)/dx;
% derivada para x,y constante
wavefrontz=gradient(phxyz)/dz+k;
%calculo de la pendiete para x,y en z=0
mzxo=wavefrontz(floor(z(1)/dz+1))/wavefrontx(N/2+1+floor(sigmaeHx(1)/dx));
mzyo=wavefrontz(floor(z(1)/dz+1))/wavefronty(N/2+1+floor(sigmaeHy(1)/dx)); 
myxo=wavefronty(N/2+1+floor(sigmaeHy(1)/dx))/wavefrontx(N/2+1+floor(sigmaeHx(1)/dx));
%-------------------Fin de calculo de Pendiente en z=0--------------------%
%-------------------------------------------------------------------------%
% Vectores para los rayos dos por cada punto en la cintura de cada uno de
% los ejes
% rayos para x
rzx1=zeros(1,sizezp(2));    
% rayos para y
rzy1=zeros(1,sizezp(2));    
% Primer posición del rayo para x
rzx1(1)=sigmaeHx(1);  
% Primer posición del rayo para y
rzy1(1)=sigmaeHy(1);

% Posición de la sombra en la obstrucción en z=0
rzx1=sigmaeHx(1);   
rzy1=sigmaeHy(1);   

% Pendientes calculadas en z=0
mzx=mzxo;
mzy=mzyo;

rx1=sigmaeHx(1);
ry1=sigmaeHy(1);

for ii=2:sizezp(2) %corriendo todos los valores de sp
    %------------Calculo de rayos-----------------------------------------%  
    % Dado z encontramos z+dz y la posición de r que le corresponde a este
    % punto z+dz con la pendiente del frente de onda que hay entre estos 
    % dos puntos
    % para x
    rx1=(1/mzx)*(z(ii)-z(ii-1))+rx1;
    % para y
    ry1=(1/mzy)*(z(ii)-z(ii-1))+ry1;
    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta
    % para x
    rzx1(ii)=rx1;   
    % para y
    rzy1(ii)=ry1;   
    % Las dos soluciones del Elegant Hermite para x e y para z=cte
    xi=x;
    yi=x;
    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=C*(qo./qz).^((nu+1/2));
    Am=C*(qo./qz).^((mu+1/2));
    An=Am;
    alphax=sqrt(i*(k)./(2*qz)).*xi;
    alphay=sqrt(i*(k)./(2*qz)).*yi;
    [eHGzx,eNHGzx]=ehermite(nu,(alphax));
    [eHGzy,eNHGzy]=ehermite(mu,(alphay));
    eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
    eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy;
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
    %-------------------Pendiente en x=sigmaeHx(1),y=0,z=ii--------------------%
    %-----------------------eHermite para x,y,z-------------------------------%
    %-----------Elegant Hermite para x,z=cte y todo y-------------------------%
    xi=rx1;     yi=x;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+1i*zl;    C=1;
    An=C*(qo./qz).^((nu+1/2));
    Am=C*(qo./qz).^((mu+1/2));
    alphax=sqrt(1i*(k)./(2*qz)).*xi;
    alphay=sqrt(1i*(k)./(2*qz)).*yi;
    [eHGxzx,eNHGxzx]=ehermite(nu,(alphax));
    [eHGxzy,eNHGxzy]=ehermite(mu,(alphay));
    eHGxzx=An.*eHGzx; eNHGxzx=An.*eNHGxzx;
    eHGxzy=Am.*eHGzy; eNHGxzy=Am.*eNHGxzy; 
    %Funciones Hankel de elegant Hermite 1D para x,z=cte y todo y%
    Hxz1x=eHGxzx+1i*eNHGxzx;
    Hxz2x=eHGxzx-1i*eNHGxzx;
    Hxz1y=eHGxzy+1i*eNHGxzy;
    Hxz2y=eHGxzy-1i*eNHGxzy;
    %Funciones Hankel de elegant Hermite 2D para x,z=cte%
    Hxz1yx=(eHGxzy).*Hxz1x;
    Hxz2yx=(eHGxzy).*Hxz2x;
    Hxz1xy= (Hxz1y).*eHGxzx;
    Hxz2xy= (Hxz2y).*eHGxzx;
    %-----------Elegant Hermite para y,z=cte y todo x-------------------------%
    xi=x;     yi=ry1;    zi=z(ii);
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=C*(qo./qz).^((nu+1/2));
    Am=C*(qo./qz).^((mu+1/2));
    alphax=sqrt(i*(k)./(2*qz)).*xi;
    alphay=sqrt(i*(k)./(2*qz)).*yi;
    [eHGyzx,eNHGyzx]=ehermite(nu,(alphax));
    [eHGyzy,eNHGyzy]=ehermite(mu,(alphay));
    eHGyzx=An.*eHGyzx; eNHGyzx=An.*eNHGzx;
    eHGyzy=Am.*eHGyzy; eNHGyzy=Am.*eNHGzy; 
    %Funciones Hankel de elegant Hermite 1D para y,z=cte%
    Hyz1x=eHGyzx+i*eNHGyzx;
    Hyz2x=eHGyzx-i*eNHGyzx;
    Hyz1y=eHGyzy+i*eNHGyzy;
    Hyz2y=eHGyzy-i*eNHGyzy;
    %Funciones Hankel de elegant Hermite 2D para y,z=cte%
    Hyz1yx=(eHGyzy).*Hyz1x;
    Hyz2yx=(eHGyzy).*Hyz2x;
    Hyz1xy=(Hyz1y).*eHGyzx;
    Hyz2xy=(Hyz2y).*eHGyzx;
    %-----------Elegant Hermite para x,y=cte y todo z-------------------------%
    xi=rx1; yi=ry1;   zi=z;
    % Parametros de eHermite
    qz=zi+i*zl; 
    An=C*(qo./qz).^((nu+1/2));
    Am=C*(qo./qz).^((mu+1/2));
    alphax=sqrt(i*(k)./(2*qz)).*xi;
    alphay=sqrt(i*(k)./(2*qz)).*yi;
    [eHGxz,eNHGxz]=ehermite(nu,(alphax));
    [eHGyz,eNHGyz]=ehermite(mu,(alphay));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    %Funciones Hankel de elegant Hermite 1D x,y=cte para todo z%
    Hx1z=eHGxz+i*eNHGxz;
    Hx2z=eHGxz-i*eNHGxz;
    Hy1z=eHGyz+i*eNHGyz;
    Hy2z=eHGyz-i*eNHGyz;
    %Funciones Hankel de elegant Hermite 2D x,y=cte para todo z%
    Hyx1z=(eHGyz).*Hx1z;
    Hyx2z=(eHGyz).*Hx2z;
    Hxy1z=(Hy1z).*eHGxz;
    Hxy2z=(Hy2z).*eHGxz;
    % Calculando la pendiente con las derivadas de las fases
    % fase x,z constante y todo y
    phxzy=unwrap(angle(Hxz2xy));
    % fase y,z constante y todo x
    phyzx=unwrap(angle(Hyz2yx));
    % fase x,y constante y todo z
    phxyz=unwrap(angle(Hxy2z));
    % derivada para x,z constante 
    wavefronty=gradient(phxzy)/dx;
    % derivada para y,z constante 
    wavefrontx=gradient(phyzx)/dx;
    % derivada para x,y constante
    wavefrontz=gradient(phxyz)/dz+k;
    %calculo de la pendiete para x,y en z=0
    mzx=wavefrontz(floor(z(ii)/dz+1))/wavefrontx(N/2+1+floor(rx1/dx));
    mzy=wavefrontz(floor(z(ii)/dz+1))/wavefronty(N/2+1+floor(ry1/dx)); 
    myx=wavefronty(N/2+1+floor(ry1/dx))/wavefrontx(N/2+1+floor(rx1/dx));
    %------------------termina calculo de los rayos-----------------------%
    figure(10)
    fig10=figure(10);
    fig10.Name=([' z = ',num2str(z(ii)/(sqrt(2)*wo))]);
    imagesc(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs((eHGzy')*(eHGzx)))
    axis square
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis1=gca;
    set(axis1,'FontSize',13);
    xlabel('$x$','Interpreter','latex','FontSize',18) 
    ylabel('$y$','Interpreter','latex','FontSize',18)
    hold on
    %lineas de la cintura en x,y
    line([ sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
    line([-sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
    line([ sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [ sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
    line([-sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo),-sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineStyle','--','LineWidth',1.5)
    %dibujando rectas de los frentes de onda
    % Ondas viajeras out
    line([ ((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo), ((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([-((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo),-((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([ ((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo),-((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [ ((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([-((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo), ((1/mzxo)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo),-((1/mzyo)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    % Ondas viajeras in-out
    line([ ((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo), ((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([-((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo),-((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([ ((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo),-((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [ ((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([-((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo), ((1/mzxo)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo),-((1/mzyo)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    % Dibujando circulo de valor esperado
    re=vesperado(ii)/(sqrt(2)*wo);
    xp=re*cos(ang);
    yp=re*sin(ang);
    plot(xp,yp,'LineWidth',1.5,'color','m');
    % Valores esperados en x, y
    line([ vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([-vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([ vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [ vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([-vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo),-vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    plot(rx1/(sqrt(2)*wo),ry1/(sqrt(2)*wo),'*','color','r')
    hold off
    pause(.001)
end  
    

