function [mzr]=gradientrz(fr,fz,k,dx,dz,x,z) 
%funcion que calcula el gradiente en un punto dado 
% f(x,y,z) a y,z constante en la primer componente
% f(x,y,z) a x,z constante en la segunda componente
% f(x,y,z) a x,y constante en la tercer componente

%derivadas parciales
gz=gradient(fr)/dx;
gr=gradient(fz)/dz+k;
N=size(gr,2);
%pendientes
mzr=gz(floor(z/dz)+1)/gr(N/2+1+floor(x/dx));



