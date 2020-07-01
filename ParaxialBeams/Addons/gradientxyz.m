function [rayObject]=gradientxyz(fx,fy,fz,k,difr,rayObject) 
%% Function for estimate gradient in (x,y,z) point of ray 
% fx = f(x,y,z) to y,z constants
% fy = f(x,y,z) to x,z constants
% fz = f(x,y,z) to x,y constants

  % point to evaluate gradient
  x =  rayObject.xCoordinate;
  y =  rayObject.yCoordinate;
  z =  rayObject.zCoordinate;
  
  % each component of diferential dr 
  tempdr     = num2cell(difr);
  [dx,dy,dz] = deal(tempdr{:});

  % partial derivatives
  gx=gradient(fx)/dx;
  gy=gradient(fy)/dy;
  gz=gradient(fz)/dz+k;
  N=size(gx,2);
  
  % get slopes 
  rayObject.zxSlope = gz(floor(z/dz+1))    /gx(N/2+1+floor(x/dx));
  rayObject.zySlope = gz(floor(z/dz+1))    /gy(N/2+1+floor(y/dy));
  rayObject.zySlope = gx(N/2+1+floor(x/dx))/gy(N/2+1+floor(y/dy));
  

end

