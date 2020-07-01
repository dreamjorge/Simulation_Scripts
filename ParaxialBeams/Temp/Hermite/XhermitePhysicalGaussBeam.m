function NHG = XhermitePhysicalGaussBeam(nu,wo,zo,x,z)

  wz   = waistPhysicalGaussianBeam(z,wo,zo);

[~,NH] = hermiteSolutions(nu,sqrt(2)*x./wz);

NHG = NH.*physicalRadialGaussianBeam(wo,zo,x,z);

end