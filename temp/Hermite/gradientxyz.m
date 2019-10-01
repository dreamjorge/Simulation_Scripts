function [mzx,mzy,mxy]=gradientxyz(fyz,fxz,fxy,k,dx,dy,dz,x,y,z) 
%funcion que calcula el gradiente en un punto dado 
% f(x,y,z) a y,z constante en la primer componente
% f(x,y,z) a x,z constante en la segunda componente
% f(x,y,z) a x,y constante en la tercer componente

%derivadas parciales
gx=gradient(fyz)/dx;
gy=gradient(fxz)/dy;
gz=gradient(fxy)/dz+k;
Ng=size(gx);
N=Ng(2);
%pendientes
mzx=gz(floor(z/dz+1))/gx(N/2+1+floor(x/dx));
mzy=gz(floor(z/dz+1))/gy(N/2+1+floor(y/dy));
mxy=gx(N/2+1+floor(x/dx))/gy(N/2+1+floor(y/dy));



