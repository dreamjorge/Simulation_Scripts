%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=28
% Parametros fisicos [micras]
wo=0.02644;
lamb=0.00006328;
k=2*pi/lamb;
zl=k*wo^2/2;
qo=-i*zl;
% El zls(elegante) es zl/2
%------------------------muestreo de vectores-----------------------------%
N=2^10;                     % Número de puntos para x, y
n=-N/2+0.05:N/2-1+0.05;     % Vector indicial igualmente espaciado

% Muestreo de vector z
Dz=zl;%/5; % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^8;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=0:dz:Dz;       % Vector z de propagacion

% Cintura de gaussiana a la distancia del tamaño de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2);
% Cintura de Laguerre a la distancia del tamaño de la ventana en z
niu=nu;
sigmaH=ws*sqrt((2*niu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tamaño de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(sigmaH)*1.8;%*1.6/(1.97)^2;  % Tamaño de la ventana
dx=Dx/N;            % Resolución
x=n*dx;             % Vector x (plano de campo a propapar)
y=x;
% Muestreo del vector de frecuencia
Du=1/dx;        % Tamaño de la ventana de Fourier
du=1/Dx;        % Resolución
u=n*du;         % Vector
% Vectores kx,ky
kx=2*pi*u;

%-----------------Elegant Hermite sin obstrucción en z=0------------------%

%Las dos soluciones del Elegant Hermite para x e y

alpha=sqrt(k*x.^2/zl);

[HGx,NHGx]=ehermite(nu,(alpha));

%funciones Hankel 1 en x
H1x=HGx+i*NHGx;

% Función a propagar
g=HGx;

%
%-------------------------Propagagación física----------------------------%

%propagador paraxial
prop=exp(-i*pi*lamb*dz*(u.^2));
figure(3)
plot(u,(angle(prop)))
title(['Propagador'])

% calculando el tamaño de zp
sizezp=size(z);

%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

%calculo de la cintura en z=0
sigmaeH=zeros(1,sizezp(2));     %vector de zeros
Iy=g.*conj(g);                   %intensidad
sigmaeH(1)=sqrt(2*(trapz(x,Iy'.*(x.^2)))./(trapz(x,Iy')));
%ciclo de propagacion
for ii=2:sizezp(2) %corriendo todos los valores de zp a excepcion de z=0
    
    %Transformada de fourier de campo a propagar
    G=fftshift(fft(g));
    %obteniendo el campo propagado
    g=ifft(fftshift(G.*prop));

    %campo g propagado
    fig4=figure(4);
    fig4.Name=([' z = ',num2str(z(ii)/(sqrt(2)*wo))]);
    pxy=g(N,N); g(N,N)=0.0845; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
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
    fig5=figure(5);
%   Para HL2 (para que estén normalizados HL1 y Hl2)
    pxy=gy(1,1); gy(1,1)=0.0845; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
    imagesc(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
    fig5.Name=([' z = ',num2str(z(ii)/(zl/2))]);
    colormap(mapgreen)
    set(gca,'YDir','normal')
    gy(1,1)=pxy;
    brighten(0.45)
    pbaspect([2.5 1 2])
    axis1=gca;
    set(axis1,'FontSize',13);
    xlabel('$z$','Interpreter','latex','FontSize',18) 
    ylabel('$x$','Interpreter','latex','FontSize',18) 
    
    %Calculo de cintura "sigmaeH"
    Iy=g(:,N/2+1).*conj(g(:,N/2+1));                   %intensidad
    sigmaeH(ii)=sqrt(2*(trapz(x,Iy'.*(x.^2)))./(trapz(x,Iy')));
    %graficando campo lateral junto con cintura
    figure(6)
    plot(x/(sqrt(2)*wo),abs(g(:,N/2+1)))
    line([sigmaeH(ii)/(sqrt(2)*wo),sigmaeH(ii)/(sqrt(2)*wo)], ylim,'color','k');
    line([-sigmaeH(ii)/(sqrt(2)*wo),-sigmaeH(ii)/(sqrt(2)*wo)], ylim,'color','k');
    title(['Propagacion Lateral ', 'z = ',num2str(z(ii)/(zl/2))])

    pause(0.001)

end

%Calculando las pendientes a los frente de ondas en la cintura del Elegant
%Laguerre en z=0

g=angle(H2(:,N/2+1));
phasedr=diff(g')./diff(x);
phasedz=k;       % Siempre es 2 (k=2)  
m1=phasedz/phasedr(floor((sigmaeH(1))/dx)+N/2+1);  % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
m2=phasedz/phasedr(floor(-(sigmaeH(1))/dx)+N/2+1);  % Pendiente 2             

%Campo lateral junto con la curva de la cintura y los frentes e onda
%normales para z=0 en la cintura del Elegant Laguerre
figure(7)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.25)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),((1/m1)*z+sigmaeH(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z-sigmaeH(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m1)*z-sigmaeH(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z+sigmaeH(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
plot(z/(zl/2),sigmaeH/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeH/(sqrt(2)*wo),'b','LineWidth',1.5)
hold off
pbaspect([2.5 1 2])
