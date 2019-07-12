%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%
mapgreen = AdvancedColormap('kgg',256,[0 150 255]/255);  %color del haz
ro=[255,69,0]/255;
%-------------------Indices del Elegant Laguerre--------------------------%
nu=4;mu=1;
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
Dz=zl*2;           % Tamaño de la ventana en z (distancia a la cual propagar)
Nz=2^6;          % Número de puntos en z
dz=Dz/Nz;        % Resolucion en z
z=0:dz:Dz;       % Vector z de propagacion
% calculando el tamaño de zp
sizezp=size(z);

% Cintura de gaussiana a la distancia del tamaño de la ventana en z
ws=wo*sqrt(1+Dz.^2/zl^2);
% Cintura de Laguerre norma la la distancia del tamaño de la ventana en z
sigmaL=ws*sqrt((2*nu+mu+1));

% Muestreo de los vectores espaciales x,y
% Calculamos el tamaño de ventana en x,y en terminos de la cintura a la
% distancia z que se propaga (esto para que no se salga de la ventana)
Dx=(2*sigmaL)*1.8;%*0.79; % Tamaño de la ventana
dx=Dx/N;            % Resolución
x=n*dx;             % Vector x (plano de campo a propapar)
y=x;
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
% Para dibujar circulos 
ang=0:0.01:2*pi; 

%---------------Elengant Laguerre sin obstrucción en z=0------------------%

% Las dos soluciones del Elegant Laguerre
Ln=exp(i*mu*TH).*R.^(abs(mu)).*Laguerre(nu,abs(mu),-i*k*R.^2/(2*qo)).*exp(i*k*R.^2/(2*qo));
Xn=exp(i*mu*TH).*R.^(abs(mu)).*XLaguerre(77,nu,abs(mu),-i*k*R.^2/(2*qo)).*exp(i*k*R.^2/(2*qo));

% Funciones Hankel Elegant Laguerre
H1=Ln+i*Xn;
H2=Ln-i*Xn;

g=Ln;

%calculo de la cintura del Elegant Laguerre en z=0
vesperado=zeros(1,sizezp(2));
sigmaeL=zeros(1,sizezp(2));     %vector de zeros
I=g.*conj(g);                   %intensidad
sigmaeL(1)=sqrt(2*(trapz(y,trapz(x,I.*(R.^2),2)))./(trapz(y,trapz(x,I,2))));
vesperado(1)=((trapz(y,trapz(x,I.*(R),2)))./(trapz(y,trapz(x,I,2))));

%Calculando las pendientes a los frente de ondas en la cintura del Elegant
%Laguerre en z=0

gan=angle(H2(:,N/2+1));
phasedr=diff(gan')./diff(x);
phasedz=k;       % Siempre es 2 (k=2)  
m1=phasedz/phasedr(floor((sigmaeL(1))/dx)+N/2+1);  % Pendiente 1 $dr((lo+xt)/0.01)$ derivada evaluada en lo+xt
m2=phasedz/phasedr(floor(-(sigmaeL(1))/dx)+N/2+1);  % Pendiente 2             



% Graficando el Elegant Laguerre en z=0
% figure(1)
% pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g).^.75)
% axis square
% shading flat
% colormap(mapgreen)
% axis1=gca;
% set(axis1,'FontSize',21);
% xlabel('$x$','Interpreter','latex','FontSize',24) 
% ylabel('$y$','Interpreter','latex','FontSize',24) 
% xlim([-2 2])
% ylim([-2 2])
% hold on
% %dibujando circulo de cintura
% ra=sigmaeL(1)/(sqrt(2)*wo);  %radio del circulo es el tamaño de ka cintura radial
% xp=ra*cos(ang);
% yp=ra*sin(ang);
% plot(xp,yp,'LineWidth',1.5,'color','r');
% %dibujando circulo de valor esperado
% re=vesperado(1)/(sqrt(2)*wo);
% xp=re*cos(ang);
% yp=re*sin(ang);
% plot(xp,yp,'LineWidth',1.5,'color','c');
% %dibujando circulo de frentes de onda
% rf1=((1/m1)*z(1)+sigmaeL(1))/(sqrt(2)*wo);
% xp=rf1*cos(ang);
% yp=rf1*sin(ang);
% plot(xp,yp,'LineWidth',1.5,'color','y');
% title(['Elegant Laguerre z=0'])
% hold off
% %Grafica lateral
% figure(2)
% plot(x/(sqrt(2)*wo),abs(g(N/2+1,:)))
% title(['Laguerre z=0 lateral'])

% %---------------Funcion a propagar con obstrucción en z=0-----------------%
% lo=0.005;               %tamaño de la obstrucción (radio de la obstruccion)
% xt=0; yt=0.0025; %xt=0.75 %traslado de la obstruccion en coordenadas cartesianas
% [thetao,ro]=cart2pol(X-xt,X'-yt);     %aplicando la traslación a coordendas xy
% obo=double(ro<=lo);  %creando la obstrucción   
% clear thetao ro      %limpienado coordenadas trasladadas
% % Función a propagar con obstrucción
% g=g.*(1-obo);
% % Recuerde normalizar respecto al máximo de la función en cuestión
% pxy=g(1,1); g(1,1)=1.9189; % Para H2 
% % Graficando función a propagar con obstrucción
% figure(2)
% pcolor(x/(sqrt(2)*wo),x/(sqrt(2)*wo),abs(g))
% axis square
% shading flat
% colormap(mapgreen)
% axis1=gca;
% set(axis1,'FontSize',13);
% xlabel('$x$','Interpreter','latex','FontSize',18)
% ylabel('$y$','Interpreter','latex','FontSize',18)
% brighten(.45)
% g(1,1)=pxy;

%-------------------------Propagagación física----------------------------%

%propagador paraxial
prop=exp(-i*pi*lamb*dz*(U.^2+(U').^2));
figure(3)
imagesc(u,u,(angle(prop)))
title(['Propagador'])



%matrices de campos transversales para guardar los datos
gx=zeros(N,sizezp(2)); gy=zeros(N,sizezp(2));

%guardando campo transversal en z=0
gx(:,1)=g(N/2+1,:);
gy(:,1)=g(:,N/2+1);


%hacer pelicula
vidObj1 = VideoWriter('circulos.avi');
vidObj1.Quality = 100;
vidObj1.FrameRate = 10;
open(vidObj1);
%ciclo de propagacion
for ii=2:sizezp(2) %corriendo todos los valores de zp a excepcion de z=0
    
    %Transformada de fourier de campo a propagar
    G=fftshift(fft2(g));
    %obteniendo el campo propagado
    g=ifft2(fftshift(G.*prop));
    %cuardando campo transversal
    gx(:,ii)=g(N/2+1,:);
    gy(:,ii)=g(:,N/2+1);

    %calculo de cintura de Elegant Laguerre en z=zk
    I=g.*conj(g);
    sigmaeL(ii)=sqrt(2*(trapz(y,trapz(x,I.*(R.^2),2)))./(trapz(y,trapz(x,I,2))));
    vesperado(ii)=((trapz(y,trapz(x,I.*(R),2)))./(trapz(y,trapz(x,I,2))));
    %campo g propagado
    fig4=figure(4);
    fig4.Name=([' z = ',num2str(z(ii)/(zl/2))]);
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
%     %dibujando circulo de la cintura
%     ra=sigmaeL(ii)/(sqrt(2)*wo);  %radio del circulo es el tamaño de ka cintura radial
%     xp=ra*cos(ang);
%     yp=ra*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','r');
    %dibujando circulo de valor esperado
%     re=vesperado(ii)/(sqrt(2)*wo); 
%     xp=re*cos(ang);
%     yp=re*sin(ang);
%     plot(xp,yp,'LineWidth',1.5,'color','y');
    %dibujando circulos de frentes de onda
    rf1=((1/m2)*z(ii)+sigmaeL(1))/(sqrt(2)*wo);
    xp=rf1*cos(ang);
    yp=rf1*sin(ang);
    plot(xp,yp,'LineWidth',2.5,'color','c');
    rf2=(-(1/m2)*z(ii)+sigmaeL(1))/(sqrt(2)*wo);
    xp=rf2*cos(ang);
    yp=rf2*sin(ang);
    plot(xp,yp,'LineWidth',2.5,'color',ro);
    hold off
    writeVideo(vidObj1, getframe(gca));
    
    %campo g en y,z
%     fig5=figure(5);
% %   Para HL2 (para que estén normalizados HL1 y Hl2)
%     pxy=gy(1,1); gy(1,1)=0.0845; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
%     imagesc(z/(zl/2),x/(sqrt(2)*wo),abs(gy))
%     fig5.Name=([' z = ',num2str(z(ii)/(zl/2))]);
%     colormap(mapgreen)
%     set(gca,'YDir','normal')
%     gy(1,1)=pxy;
%     brighten(0.45)
%     pbaspect([2 1 2])
%     axis1=gca;
%     set(axis1,'FontSize',13);
%     xlabel('$z$','Interpreter','latex','FontSize',18) 
%     ylabel('$x$','Interpreter','latex','FontSize',18) 
%     %graficando las rectas de las normales al frente de onda
%     hold on
%     plot(z(1:ii)/(zl/2),((1/m1)*z(1:ii)+sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
%     plot(z(1:ii)/(zl/2),((1/m2)*z(1:ii)-sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
%     plot(z(1:ii)/(zl/2),((1/m1)*z(1:ii)-sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%     plot(z(1:ii)/(zl/2),((1/m2)*z(1:ii)+sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%     plot(z(1:ii)/(zl/2),vesperado(1:ii)/(sqrt(2)*wo),'y','LineWidth',1.5)
%     plot(z(1:ii)/(zl/2),-vesperado(1:ii)/(sqrt(2)*wo),'y','LineWidth',1.5)
% 	writeVideo(vidObj2, getframe(gca));
%     hold off
   
%     %graficando Elegant Laguerre en lateral junto con cintura
%     figure(6)
%     plot(x/(sqrt(2)*wo),abs(g(:,N/2+1)))
%     line([sigmaeL(ii)/(sqrt(2)*wo),sigmaeL(ii)/(sqrt(2)*wo)], ylim,'color','r');
%     line([-sigmaeL(ii)/(sqrt(2)*wo),-sigmaeL(ii)/(sqrt(2)*wo)], ylim,'color','r');
%     line([vesperado(ii)/(sqrt(2)*wo),vesperado(ii)/(sqrt(2)*wo)], ylim,'color','k');
%     line([-vesperado(ii)/(sqrt(2)*wo),-vesperado(ii)/(sqrt(2)*wo)], ylim,'color','k');
%     title(['Propagacion Lateral ', 'z = ',num2str(z(ii)/(zl/2))])
    pause(0.001)

end
%# save as AVI file
close(vidObj1);
%Campo lateral junto con la curva de la cintura y los frentes e onda
%normales para z=0 en la cintura del Elegant Laguerre
figure(7)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.5)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
plot(z/(zl/2),((1/m1)*z+sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z-sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
plot(z/(zl/2),((1/m1)*z-sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
plot(z/(zl/2),((1/m2)*z+sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
% plot(z/(zl/2),sigmaeL/(sqrt(2)*wo),'b','LineWidth',1.5)
% plot(z/(zl/2),-sigmaeL/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),vesperado/(sqrt(2)*wo),'y','LineWidth',1.5)
plot(z/(zl/2),-vesperado/(sqrt(2)*wo),'y','LineWidth',1.5)
hold off
pbaspect([2.5 1 2])


figure(8)
pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gy).^.5)
shading interp
colormap(mapgreen)
% title('Campo propagado corte en x')
hold on
%graficando las rectas de las normales al frente de onda
% plot(z/(zl/2),((1/m1)*z+sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
% plot(z/(zl/2),((1/m2)*z-sigmaeL(1))/(sqrt(2)*wo),'r','LineWidth',1.5)
% plot(z/(zl/2),((1/m1)*z-sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
% plot(z/(zl/2),((1/m2)*z+sigmaeL(1))/(sqrt(2)*wo),'c','LineWidth',1.5)
%graficando la cintura del Elegant Laguerre
plot(z/(zl/2),sigmaeL/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),-sigmaeL/(sqrt(2)*wo),'b','LineWidth',1.5)
plot(z/(zl/2),vesperado/(sqrt(2)*wo),'y','LineWidth',1.5)
plot(z/(zl/2),-vesperado/(sqrt(2)*wo),'y','LineWidth',1.5)
hold off
pbaspect([2 1 2])

