%eHermite analiticos
%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=20;mu=20;
% Parametros fisicos [micras]
wo=2;
lamb=0.6328;    %micras
k=2*pi/lamb;
zl=k*wo^2/2;
qo=i*zl;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl;%/5; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^7;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=0:dz:Dz;       % Vector z de propagacion
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
%--------curvas con z=constante=0
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=x/wo;
yi=x/wo;
zi=0;
% Parametros de eHermite
qz=zi+i*zl;    C=1;
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
An=1;Am=1;
alphax=sqrt(i*(k)./(2*qz)).*xi;
alphay=sqrt(i*(k)./(2*qz)).*yi;
alphax=xi;
alphay=xi;
[eHGzx,eNHGzx]=ehermite(nu,(alphax));
[eHGzy,eNHGzy]=ehermite(mu,(alphay));
eHGzx=An.*eHGzx; eNHGzx=An.*eNHGzx;
eHGzy=Am.*eHGzy; eNHGzy=Am.*eNHGzy; 
%Funciones Hankel de elegant Hermite%
%Funcion Hankel con z=cte
H1z=(eHGzy')*(eHGzx)+i*(eNHGzy')*(eNHGzx);
H2z=(eHGzy')*(eHGzx)-i*(eNHGzy')*(eNHGzx);
figure(1)
imagesc(abs(H1z+H2z))
figure(2)
imagesc(abs(H1z+H2z))
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
%-----------------Elegant Hermite sin obstrucción en z=0------------------%
%-------------Pendientes a los frente de ondas en la cintura--------------%
%----------------------del Elegant Hermite en z=0-------------------------%
%--------------------funciones analiticas de eH para z=0------------------%
%-------------------------------------------------------------------------%
%-------curva para x,y=constante y todo z
% Las dos soluciones del Elegant Hermite para x e y para z=0
xi=sigmaeHx(1);
yi=sigmaeHy(1);
zi=z;
% Parametros de eHermite
qz=zi+i*zl; 
An=C*(qo./qz).^((nu+1/2));
Am=C*(qo./qz).^((mu+1/2));
alphax=sqrt(i*(k)./(2*qz)).*xi;
alphay=sqrt(i*(k)./(2*qz)).*xi;
[eHGxz,eNHGxz]=ehermite(nu,(alphax));
[eHGyz,eNHGyz]=ehermite(mu,(alphay));
eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
% Funciones Hankel de elegant Hermite%
% Funciones Hankel 1 en x
H1xz=eHGxz+i*eNHGxz;
H1yz=eHGyz+i*eNHGyz;
% Funciones Hankel 2 en y
H2xz=eHGxz-i*eNHGxz;
H2yz=eHGxz-i*eNHGxz;
% Funcion Hankel con x,y=cte y todo z
H1xy=(eHGyz')*(eHGxz)+i*(eNHGyz')*(eNHGxz);
H2xy=(eHGyz')*(eHGxz)-i*(eNHGyz')*(eNHGxz);

% Calculando la pendiente con las derivadas de las fases
% fase z constante y todo x,y
phzx=unwrap(angle(H2z(:,N/2+1)));
phzy=unwrap(angle(H2z(N/2+1,:)));
% % fase x,y constante y todo z
phxz=unwrap(angle(H2xy(:,Nz/2+1)));
phyz=unwrap(angle(H2xy(Nz/2+1,:)));
% fase z constante y todo x,y
% phzx=unwrap(angle(H2zx));
% phzy=unwrap(angle(H2zy));
% fase x,y constante y todo z
% phxz=unwrap(angle(H2xz));
% phyz=unwrap(angle(H2yz));
%derivadas para x,y constantes
wavefrontxz=gradient(phxz)/dz+k;
wavefrontyz=gradient(phyz)/dz+k;
%derivadas para z constante
wavefrontzx=gradient(phzx)/dx;
wavefrontzy=gradient(phzy)/dx;
%calculo de la pendiete para x,y en z=0
mxoo=wavefrontxz(floor(z(1)/dz+1))/wavefrontzx(N/2+1+floor(xi/dx));
myoo=wavefrontyz(floor(z(1)/dz+1))/wavefrontzy(N/2+1+floor(yi/dx)); 
%-------------------------------------------------------------------------%


% %--------------pendiente para z=0
% % fase z constante y todo x,y
% phzxo=unwrap(angle(H2x));
% phzyo=unwrap(angle(H2y));
% %derivadas para x,y constantes
% wavefrontxz=k;
% wavefrontyz=k;
% %derivadas para z constante
% wavefrontzxo=gradient(phzxo)/dx;
% wavefrontzyo=gradient(phzyo)/dx;
% %calculo de la pendiete para x,y en z=0
% mxoo=k/wavefrontzxo(N/2+1+floor(sigmaeHx(1)/dx));
% myoo=k/wavefrontzyo(N/2+1+floor(sigmaeHy(1)/dx)); 

%-------------------------------------------------------------------------%
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


% Vectores para los rayos dos por cada punto en la cintura de cada uno de
% los ejes
% rayos para x
rxz1=zeros(1,sizezp(2));    rxz2=zeros(1,sizezp(2));
% rayos para y
ryz1=zeros(1,sizezp(2));    ryz2=zeros(1,sizezp(2));
% Primer posición del rayo para x
rxz1(1)=sigmaeHx(1);   rxz2(1)=sigmaeHx(1);
% Primer posición del rayo para y
ryz1(1)=sigmaeHy(1);   ryz2(1)=sigmaeHy(1);

% Posición de la sombra en la obstrucción en z=0
rx1=sigmaeHx(1);   rx2=sigmaeHx(1);
ry1=sigmaeHy(1);   ry2=sigmaeHy(1);

% Pendientes calculadas en z=0
mx1=mxoo;
my1=myoo;

g=(eHGzy')*eHGzx;
% Matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));
% Guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);    gy(:,1)=g(:,N/2+1);

for ii=2:sizezp(2) %corriendo todos los valores de sp
    %------------Calculo de rayos-----------------------------------------%  
    % Dado z encontramos z+dz y la posición de r que le corresponde a este
    % punto z+dz con la pendiente del frente de onda que hay entre estos 
    % dos puntos
    % para x
    rx1=( 1/mx1)*(z(ii)-z(ii-1))+rx1;
    rx2=(-1/mx1)*(z(ii)-z(ii-1))+rx2;
    % para y
    ry1=( 1/my1)*(z(ii)-z(ii-1))+ry1;
    ry2=(-1/my1)*(z(ii)-z(ii-1))+ry2;
    % Guardando la posición de este rayo para z(ii) para generar la curva
    % a traves de segmentos de recta
    % para x
    rxz1(ii)=rx1;   rxz2(ii)=rx2;
    % para y
    ryz1(ii)=ry1;   ryz2(ii)=ry2;
    %-----------Calculando las pendientes del siguiente paso--------------%
    %--------curvas con z=constante=0
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
    %Funciones Hankel de elegant Hermite%
    % Funciones Hankel 2 en x,y con z=cte
    H2z=(eHGzy')*(eHGzx)-i*(eNHGzy')*(eNHGzy);
    %---------------------------------------------------------------------%
    %-------curva para x,y=constante y todo z
    % Las dos soluciones del Elegant Hermite para x e y para z=0
    xi=rx1;
    yi=ry1;
    zi=z;
    % Parametros de eHermite
    qz=zi+i*zl;    C=1;
    An=C*(qo./qz).^((nu+1/2));
    Am=C*(qo./qz).^((mu+1/2));
    An=1;
    An=Am;
    alphax=sqrt(i*(k)./(2*qz)).*xi;
    alphay=sqrt(i*(k)./(2*qz)).*yi;
    [eHGxz,eNHGxz]=ehermite(nu,(alphax));
    [eHGyz,eNHGyz]=ehermite(mu,(alphay));
    eHGxz=An.*eHGxz; eNHGxz=An.*eNHGxz;
    eHGyz=Am.*eHGyz; eNHGyz=Am.*eNHGyz;
    % Funciones Hankel 2 x,y=cte y todo z
    H2xy=(eHGyz')*(eHGxz)-i*(eNHGyz')*(eNHGxz);
    % ---Fases
    % ---Fases
%     % fase z constante y todo x,y
    phzx=unwrap(angle(H2z(:,N/2+1)));
    phzy=unwrap(angle(H2z(N/2+1,:)));
%     figure(10)
%     plot(phzx)
%     figure(11)
%     plot(phzy)
%     % fase x,y constante y todo z
     phxz=unwrap(angle(H2xy(:,Nz/2+1)));
     phyz=unwrap(angle(H2xy(Nz/2+1,:)));%
%     figure(12)
%     plot(phxz)
%     figure(13)
%     plot(phyz)
    % fase z constante y todo x,y
%     phzx=unwrap(angle(H2zx));
%     phzy=unwrap(angle(H2zy));
    % fase x,y constante y todo z
%     phxz=unwrap(angle(H2xz));
%     phyz=unwrap(angle(H2yz));
    % derivadas para x,y constantes
    wavefrontxz=gradient(phxz)/dz+k;
    wavefrontyz=gradient(phyz)/dz+k;
    % derivadas para z constantes
    wavefrontzx=gradient(phzx)/dx;
    wavefrontzy=gradient(phzy)/dx;

    % Pendientes
    mx1=wavefrontxz(floor(z(ii)/dz+1))/wavefrontzx(N/2+1+floor(rx1/dx));
    my1=wavefrontyz(floor(z(ii)/dz+1))/wavefrontzy(N/2+1+floor(ry1/dx));
    pause(0.001)
    g=(eHGzy')*eHGzx;
    gx(:,ii)=g(:,N/2+1);
    gy(:,ii)=g(N/2+1,:);
     %---------------Calculando las cinturas en x,y------------------------%
    Iy=g(:,N/2+1).*conj(g(:,N/2+1));	%intensidad
    Ix=g(N/2+1,:).*conj(g(N/2+1,:));	%intensidad
    I=Iy*Ix;
    sigmaeHx(ii)=sqrt(2*(trapz(x,Ix.*(x.^2)))./(trapz(x,Ix)));
    sigmaeHy(ii)=sqrt(2*(trapz(x,Iy'.*(x.^2)))./(trapz(x,Iy')));
    %------------------Calculando el valor esperado-----------------------%
    %valor esperado radial
    vesperado(ii)=((trapz(y,trapz(x,I.*sqrt(X.^2+(X').^2),2)))./(trapz(y,trapz(x,I,2))));
    % Intensidad en un cuadrante par calcular el valor esperado en las
    % componentes
    Ix1=Ix(N/2:N);  Iy1=Iy(N/2:N);  %I1=Iy1*Ix1;
    vesperadoy(ii)=(trapz(x1,Ix1.*x1))./(trapz(x1,Ix1));
    vesperadox(ii)=(trapz(x1,(Iy1').*x1))./(trapz(x1,Iy1'));
    %------------------termina calculo de los rayos-----------------------%
end  
    

figure(8)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gx).^.15)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en y')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),(-(1/myoo)*z+sigmaeHy(1))/(sqrt(2)*wo),'c--','LineWidth',1.5)
% plot(z/(zl/2),(-(1/my)*z-sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/myoo)*z-sigmaeHy(1))/(sqrt(2)*wo),'c--','LineWidth',1.5)
% plot(z/(zl/2),((1/my)*z+sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),(ryz1)/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),(ryz2)/(sqrt(2)*wo),'b','LineWidth',1.5)
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
plot(z/(zl/2),(-(1/mxoo)*z+sigmaeHx(1))/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),(-(1/mxoo)*z-sigmaeHx(1))/(sqrt(2)*wo),'r--','LineWidth',1.5)
plot(z/(zl/2),( (1/mxoo)*z-sigmaeHx(1))/(sqrt(2)*wo),'c--','LineWidth',1.5)
plot(z/(zl/2),( (1/mxoo)*z+sigmaeHx(1))/(sqrt(2)*wo),'r--','LineWidth',1.5)
plot(z/(zl/2), (ryz1)/(sqrt(2)*wo),'k','LineWidth',1.5)
plot(z/(zl/2), (ryz2)/(sqrt(2)*wo),'k','LineWidth',1.5)
plot(z/(zl/2), (-ryz1)/(sqrt(2)*wo),'k','LineWidth',1.5)
plot(z/(zl/2), (-ryz2)/(sqrt(2)*wo),'k','LineWidth',1.5)
plot(z/(zl/2), vesperadoy/(sqrt(2)*wo),'m','LineWidth',1.5)
plot(z/(zl/2),-vesperadoy/(sqrt(2)*wo),'m','LineWidth',1.5)
plot(z/(zl/2), sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
hold off
pbaspect([2.5 1 2])