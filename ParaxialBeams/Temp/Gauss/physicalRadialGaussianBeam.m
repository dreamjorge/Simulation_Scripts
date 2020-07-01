function [GB]=physicalRadialGaussianBeam(wo,zo,r,z) 
  
  k    = (2*zo)/wo^2;
  
  wz   = waistPhysicalGaussianBeam(z,wo,zo);
  Phiz = phasePhysicalGaussianBeam(z,zo);

  Rz   = radiusPhysicalGaussianBeam(z,zo);
  
  u0   = 1;
  
  GB   = (u0./wz).*exp(1i*k*(r.^2)./(2*Rz)).*exp(-1i*Phiz).*exp(-(r.^2)./((wz.^2))).*exp(1i*k*z);%%;


end