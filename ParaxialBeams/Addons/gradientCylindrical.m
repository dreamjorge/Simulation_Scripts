function [rayObject] = gradientCylindrical(fr,ftheta,fz,k,...
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
  
  grpoint = gr(floor(r/dr+1));
  gzpoint = gz(floor(z/dz+1));
  gthpoint = gtheta(N/2+floor(theta/dtheta+1));
  
  rayObject.zrSlope  = gzpoint./grpoint;
  rayObject.zthSlope = gzpoint./gthpoint;
  rayObject.rthSlope = grpoint./gthpoint;
  

end