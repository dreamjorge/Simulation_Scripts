function ws = waistPhysicalGaussianBeam(s,wo,so)

    ws = (wo)*sqrt((s./so).^2+1);

end
