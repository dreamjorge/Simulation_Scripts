function omegas = elegantWaistGaussianBeam(s,wo)

%Parameters Geometrics

so     = wo.^2;

omegas = (wo/sqrt(2))*sqrt((s./so).^2+1);

end