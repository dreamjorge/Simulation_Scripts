function [SlopesObject]=gradientxyz(fx,fy,fz,k,difr,rayObject) 
%funcion que calcula el gradiente en un punto dado 
% f(x,y,z) a y,z constante en la primer componente
% f(x,y,z) a x,z constante en la segunda componente
% f(x,y,z) a x,y constante en la tercer componente

  x          =  rayObject.xCoordinate;
  y          =  rayObject.yCoordinate;
  z          =  rayObject.zCoordinate;
  
  tempdr     = num2cell(difr);
  [dx,dy,dz] = deal(tempdr{:});

  %derivadas parciales
  gx=gradient(fx)/dx;
  gy=gradient(fy)/dy;
  gz=gradient(fz)/dz+k;
  N=size(gx,2);
  
  %pendientes
  SlopesObject    = Slopes;
  SlopesObject.zx = gz(floor(z/dz+1))    /gx(N/2+1+floor(x/dx));
  SlopesObject.zy = gz(floor(z/dz+1))    /gy(N/2+1+floor(y/dy));
  SlopesObject.xy = gx(N/2+1+floor(x/dx))/gy(N/2+1+floor(y/dy));
  

end

