function Rs = radiusGaussianBeam(s,wo)


  so = wo.^2;

  Rs = (s).*(1+(so./s).^2);

    if (find(isnan(Rs))~= 0)
        Rs(find(isnan(Rs)))= inf;
    end

end