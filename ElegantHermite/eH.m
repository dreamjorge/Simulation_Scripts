%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
%-----------------------indices del Elegant Hermite------------------------%
nu=40;mu=40;
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
Dz=zl*1.5;%/5; % Tamaño de la ventana en z (distancia a la cual propagar)
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
Dx=(sigmaH)*1.2;%*1.6/(1.97)^2;  % Tamaño de la ventana
dx=Dx/N;            % Resolución
x=n*dx;             % Vector x (plano de campo a propapar)
y=x;
[X]=meshgrid(x);
x1=x(N/2:N);        % vector de 0 a Dx
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
%-----------------Elegant Hermite sin obstrucción en z=0------------------%

% Las dos soluciones del Elegant Hermite para x e y para z=0

alpha=x/wo;

[eHGx,eNHGx]=ehermite(nu,(alpha));
[eHGy,eNHGy]=ehermite(mu,(alpha));

%-------------------Funciones Hankel de elegant Hermite-------------------%
% Funciones Hankel 1 en x
H1x=eHGx+i*eNHGx;
H1y=eHGx+i*eNHGx;
% Funciones Hankel 2 en y
H2x=eHGx-i*eNHGx;
H2y=eHGx-i*eNHGx;

% Hermite xy
eHG=(eHGy')*eHGx;
eNHG=(eNHGy')*eNHGy;

% Funcion hankel 1 xy
H1=(H1y')*H1x;
% Funcion hankel 2 xy
H2=(H2y')*H2x;


%----------------Cintura de elegant Hermite en z=0------------------------%

% Matrices de ceros para las cinturas en x,y en z=0
sigmaeHx=zeros(1,sizezp(2));        sigmaeHy=zeros(1,sizezp(2));
% Intensidades en total y en x, y en z=0
Ix=eHG(N/2+1,:).*conj(eHG(N/2+1,:)); Iy=eHG(:,N/2+1).*conj(eHG(:,N/2+1));
I=Iy*Ix;
% Calculando las cinturas
sigmaeHx(1)=sqrt(2*(trapz(x,Ix.*(x.^2)))./(trapz(x,Ix)));
sigmaeHy(1)=sqrt(2*(trapz(x,Iy'.*(x.^2)))./(trapz(x,Iy')));


%-------------Pendientes a los frente de ondas en la cintura--------------%
%----------------------del Elegant Hermite en z=0-------------------------%

ph1=angle(H2(:,N/2+1));
phasedr=diff(ph1')./diff(x);
phasedz=k;       % Siempre es 2 (k=2)  
m1=phasedz/phasedr(floor((sigmaeHx(1))/dx)+N/2+1);   % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt

ph1=angle(H2(N/2+1,:));
phasedr=diff(ph1)./diff(x);
phasedz=k;       % Siempre es 2 (k=2)  
m2=phasedz/phasedr(floor((sigmaeHy(1))/dx)+N/2+1); 

% m2=phasedz/phasedr(floor((sigmaeHy(1))/dx)+N/2+1);   % Pendiente 2  


%------------------Valor esperado radial y en componentes-----------------%

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
vesperadox(1)=(trapz(x1,(Iy1').*x1))./(trapz(x1,Iy1'));

%-------------------------Propagagación física----------------------------%

% Función a propagar
g=eHG;

% Graficando el campo en z=0
figure(1)
pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^.75)
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',21);
xlabel('$x$','Interpreter','latex','FontSize',24) 
ylabel('$y$','Interpreter','latex','FontSize',24) 
hold on
% Dibujando cinturas de x, y
line([ sigmaeHx(1)/(sqrt(2)*wo)  sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo)  sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
line([-sigmaeHx(1)/(sqrt(2)*wo) -sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo)  sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
line([ sigmaeHx(1)/(sqrt(2)*wo) -sigmaeHx(1)/(sqrt(2)*wo)], [ sigmaeHy(1)/(sqrt(2)*wo)  sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
line([-sigmaeHx(1)/(sqrt(2)*wo)  sigmaeHx(1)/(sqrt(2)*wo)], [-sigmaeHy(1)/(sqrt(2)*wo) -sigmaeHy(1)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
% Dibujando circulo de valor esperado
re=vesperado(1)/(sqrt(2)*wo);
xp=re*cos(ang);
yp=re*sin(ang);
plot(xp,yp,'LineWidth',1.5,'color','m');
% Dibujando valor esperado x,y
line([ vesperadoy(1)/(sqrt(2)*wo), vesperadoy(1)/(sqrt(2)*wo)], [-vesperadox(1)/(sqrt(2)*wo), vesperadox(1)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
line([-vesperadoy(1)/(sqrt(2)*wo),-vesperadoy(1)/(sqrt(2)*wo)], [-vesperadox(1)/(sqrt(2)*wo), vesperadox(1)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
line([ vesperadoy(1)/(sqrt(2)*wo),-vesperadoy(1)/(sqrt(2)*wo)], [ vesperadox(1)/(sqrt(2)*wo), vesperadox(1)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
line([-vesperadoy(1)/(sqrt(2)*wo), vesperadoy(1)/(sqrt(2)*wo)], [-vesperadox(1)/(sqrt(2)*wo),-vesperadox(1)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
xlim([-5 5])
ylim([-5 5])
hold off

%propagador paraxial
prop=exp(-i*pi*lamb*dz*(U.^2+(U').^2));
figure(3)
imagesc(u,u,(angle(prop)))
title(['Propagador'])

% Matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

% Guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);    gy(:,1)=g(:,N/2+1);

% Grabar videos
% vidObj1 = VideoWriter('cuadrados.avi');
% vidObj1.Quality = 100;
% vidObj1.FrameRate = 10;
% open(vidObj1);
% vidObj2 = VideoWriter('circulos.avi');
% vidObj2.Quality = 100;
% vidObj2.FrameRate = 10;
% open(vidObj2);

% Ciclo de propagacion
for ii=2:sizezp(2) %corriendo todos los valores de zp a excepcion de z=0
    
    %Transformada de fourier de campo a propagar
    G=fftshift(fft2(g));
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
    %cuardando campo transversal
    gx(:,ii)=g(:,N/2+1);
    gy(:,ii)=g(N/2+1,:);

    Iy=g(:,N/2+1).*conj(g(:,N/2+1));	%intensidad
    Ix=g(N/2+1,:).*conj(g(N/2+1,:));	%intensidad
    I=Iy*Ix;
    Ix1=Ix(N/2:N);
    Iy1=Iy(N/2:N);
    x1=x(N/2:N);
    sigmaeHx(ii)=sqrt(2*(trapz(x,Ix.*(x.^2)))./(trapz(x,Ix)));
    sigmaeHy(ii)=sqrt(2*(trapz(x,Iy'.*(x.^2)))./(trapz(x,Iy')));
    
    %calculando el valor esperado
    vesperado(ii)=((trapz(y,trapz(x,I.*sqrt(X.^2+(X').^2),2)))./(trapz(y,trapz(x,I,2))));
    vesperadoy(ii)=(trapz(x1,Ix1.*x1))./(trapz(x1,Ix1));
    vesperadox(ii)=(trapz(x1,(Iy1').*x1))./(trapz(x1,Iy1'));
    
    %campo g propagado con cuadrados
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
    hold on
    %lineas de la cintura en x,y
    line([ sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
    line([-sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
    line([ sigmaeHx(ii)/(sqrt(2)*wo),-sigmaeHx(ii)/(sqrt(2)*wo)], [ sigmaeHy(ii)/(sqrt(2)*wo), sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
    line([-sigmaeHx(ii)/(sqrt(2)*wo), sigmaeHx(ii)/(sqrt(2)*wo)], [-sigmaeHy(ii)/(sqrt(2)*wo),-sigmaeHy(ii)/(sqrt(2)*wo)],'color','b','LineWidth',1.5)
    %dibujando rectas de los frentes de onda
    %out
    line([ ((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo), ((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([-((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo),-((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([ ((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo),-((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [ ((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    line([-((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo), ((1/m1)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo),-((1/m2)*z(ii)-sigmaeHy(1))/(sqrt(2)*wo)],'color','r','LineWidth',1.5)
    %in-out
    line([ ((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo), ((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([-((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo),-((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([ ((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo),-((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [ ((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo), ((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
    line([-((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo), ((1/m1)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo)], [-((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo),-((1/m2)*z(ii)+sigmaeHy(1))/(sqrt(2)*wo)],'color','c','LineWidth',1.5)
        %dibujando circulo de valor esperado
    re=vesperado(ii)/(sqrt(2)*wo);
    xp=re*cos(ang);
    yp=re*sin(ang);
    plot(xp,yp,'LineWidth',1.5,'color','m');
    % Valores esperados en x, y
    line([ vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([-vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([ vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [ vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    line([-vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo),-vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
    hold off
%     writeVideo(vidObj1, getframe(gca));

%     %campo g propagado con circulos
%     fig41=figure(41);
%     fig41.Name=([' z = ',num2str(z(ii)/(sqrt(2)*wo))]);
%     pxy=g(N,N); g(N,N)=0.0845; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
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
%     hold on
%     %dibujando circulo de valor esperado
%     re=vesperado(ii)/(sqrt(2)*wo);
%     xp=re*cos(ang);
%     yp=re*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','m');
%     % Valores esperados en x, y
%     line([ vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([-vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([ vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], [ vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     line([-vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], [-vesperadox(ii)/(sqrt(2)*wo),-vesperadox(ii)/(sqrt(2)*wo)],'color','m','LineWidth',1.5)
%     %circulo out
%     rout=sqrt(2)*((1/m2)*z(ii)+sigmaeHx(1))/(sqrt(2)*wo);
%     xp=rout*cos(ang);
%     yp=rout*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','r');
%     %circulo in-out
%     rin=sqrt(2)*((1/m2)*z(ii)-sigmaeHx(1))/(sqrt(2)*wo);
%     xp=rin*cos(ang);
%     yp=rin*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','c');
%     hold off
% 	writeVideo(vidObj2, getframe(gca));


%     %campo g en y,z
%     fig5=figure(5);
% %   Para HL2 (para que estén normalizados HL1 y Hl2)
%     pxy=gy(1,1); gy(1,1)=0.0845; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
%     imagesc(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
%     fig5.Name=([' z = ',num2str(z(ii)/(zl/2))]);
%     colormap(mapgreen)
%     set(gca,'YDir','normal')
%     gy(1,1)=pxy;
%     brighten(0.45)
%     pbaspect([2.5 1 2])
%     axis1=gca;
%     set(axis1,'FontSize',13);
%     xlabel('$z$','Interpreter','latex','FontSize',18) 
%     ylabel('$x$','Interpreter','latex','FontSize',18) 
    
    
    %graficando campo lateral junto con el valor esperado
    figure(6)
    plot(x/(sqrt(2)*wo),unwrap(angle(g(:,N/2+1))))
    line([ vesperadox(ii)/(sqrt(2)*wo), vesperadox(ii)/(sqrt(2)*wo)], ylim,'color','k');
    line([-vesperadox(ii)/(sqrt(2)*wo),-vesperadox(ii)/(sqrt(2)*wo)], ylim,'color','k');
    title(['Propagacion Lateral x  ', 'z = ',num2str(z(ii)/(zl/2))])

    figure(7)
    plot(x/(sqrt(2)*wo),unwrap(angle(g(N/2+1,:))))
    line([ vesperadoy(ii)/(sqrt(2)*wo), vesperadoy(ii)/(sqrt(2)*wo)], ylim,'color','k');
    line([-vesperadoy(ii)/(sqrt(2)*wo),-vesperadoy(ii)/(sqrt(2)*wo)], ylim,'color','k');
    title(['Propagacion Lateral y  ', 'z = ',num2str(z(ii)/(zl/2))])
    
    pause(0.001)

end
%# save as AVI file
% close(vidObj1);close(vidObj2);

%Campo lateral junto con la curva de la cintura y los frentes e onda
%normales para z=0 en la cintura del Elegant Laguerre
figure(8)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gx).^.15)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),(-(1/m1)*z+sigmaeHx(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),(-(1/m1)*z-sigmaeHx(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m1)*z-sigmaeHx(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),((1/m1)*z+sigmaeHx(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
plot(z/(zl/2),sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeHx/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),vesperadox/(sqrt(2)*wo),'m','LineWidth',1.5)
plot(z/(zl/2),-vesperadox/(sqrt(2)*wo),'m','LineWidth',1.5)
hold off
pbaspect([2.5 1 2])

figure(9)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.15)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en y')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),(-(1/m2)*z+sigmaeHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),(-(1/m2)*z-sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z-sigmaeHy(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z+sigmaeHy(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
plot(z/(zl/2),sigmaeHy/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeHy/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),vesperadox/(sqrt(2)*wo),'m','LineWidth',1.5)
plot(z/(zl/2),-vesperadox/(sqrt(2)*wo),'m','LineWidth',1.5)
hold off
pbaspect([2.5 1 2])


