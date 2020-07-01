function HG = hermitePhysicalGaussBeam(nu,wo,zo,x,z)

[H,~] = hermiteSolutions(nu,x);

HG = H.*physicalRadialGaussianBeam(wo,zo,x,z);