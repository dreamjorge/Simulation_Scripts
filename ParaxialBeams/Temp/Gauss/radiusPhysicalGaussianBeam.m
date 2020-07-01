function Rs = radiusPhysicalGaussianBeam(s,so)

        Rs = (s).*(1+(so./s).^2);
        
        if (find(isnan(Rs))~= 0)
            Rs(find(isnan(Rs)))= inf;
        end
end