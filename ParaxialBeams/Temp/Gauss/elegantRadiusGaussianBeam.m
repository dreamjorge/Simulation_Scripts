function Rs = elegantRadiusGaussianBeam(s,wo)


so = wo.^2;

Rs = (s/2).*(1+(so./s).^2);

end