function ws = waistGaussianBeam(s,wo)

%Parameters Geometrics

so     = wo.^2;

ws = (wo)*sqrt((s./so).^2+1);

end
