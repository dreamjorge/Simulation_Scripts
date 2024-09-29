function HG = hermitePhysicalGaussBeam(nu,wo,zo,x,z)

  wz   = waistPhysicalGaussianBeam(z,wo,zo);

[H,~] = hermiteSolutions(nu,sqrt(2)*x./wz);

HG = H.*physicalRadialGaussianBeam(wo,zo,x,z);

end