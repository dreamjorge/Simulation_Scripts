function [SlopesObject] = gradientCylindrical(fr,ftheta,fz,k,...
                                           difr,...
                                           rayObject) 

  r      =  rayObject.rCoordinate;
  theta  =  rayObject.thetaCoordinate;
  z      =  rayObject.zCoordinate;
  
  tempdr         = num2cell(difr);
  [dr,dtheta,dz] = deal(tempdr{:});


  %partial derivatives 
  gr     = gradient(fr)/dr;
  gtheta = (1./r).*(gradient(ftheta)/dtheta);
  gz     = gradient(fz)/dz+k;
  N      = size(gr,2);

  SlopesObject = Slopes;
  SlopesObject.zr  = gz(floor(z/dz+1))    /gr(N/2+1+floor(r/dr));
  SlopesObject.zth = gz(floor(z/dz+1))    /gtheta(N/2+1+floor(theta/dtheta));
  SlopesObject.rth = gr(N/2+1+floor(r/dr))/gtheta(N/2+1+floor(theta/dtheta));
  

end