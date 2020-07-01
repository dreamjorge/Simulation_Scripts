function sigmas = sigmaGaussianBeam(s,wo)

%Parameters Geometrics

so     = wo.^2;

sigmas=(wo/2)*sqrt((s./so).^2+1);

end