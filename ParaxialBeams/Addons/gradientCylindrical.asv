function [mzr,mztheta,mrtheta] = gradientCylindrical(fr,ftheta,fz,k,dr,dtheta,dz,r,theta,z) 

  %partial derivatives
  
  gr     = gradient(fr)/dr;
  gtheta = (1./r).*(gradient(ftheta)/dtheta);
  gz     = gradient(fz)/dz+k;
  N      = size(gr,2);

  %
  mzr     = gz(floor(z/dz+1))    /gr(N/2+1+floor(r/dr));
  mztheta = gz(floor(z/dz+1))    /gtheta(N/2+1+floor(theta/dtheta));
  mrtheta = gr(N/2+1+floor(r/dr))/gtheta(N/2+1+floor(theta/dtheta));

end