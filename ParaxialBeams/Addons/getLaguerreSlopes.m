function [ray] = getLaguerreSlopes(ray,x,y,z,...
                                       dx,dy,dz,...
                                       xi,yi,zi,...
                                       InitialWaist,Wavelength,p,l,nh)

  %y = cte, z = cte
  HLx = HankelLaguerre(x,yi,zi,InitialWaist,Wavelength,p,l,nh);

  %x = cte, z = cte
  HLy = HankelLaguerre(xi,y,zi,InitialWaist,Wavelength,p,l,nh);

  %x = cte, y = cte
  HLz = HankelLaguerre(xi,yi,z,InitialWaist,Wavelength,p,l,nh);

  fx   = unwrap(angle(HLx.OpticalField));
  fy   = unwrap(angle(HLy.OpticalField));
  fz   = unwrap(angle(HLz.OpticalField));

  k = 2*pi/InitialWaist;
  
  [ray.zxSlope,ray.zySlope,ray.xySlope] = gradientxyz(fx,fy,fz,k,dx,dy,dz,xi,yi,zi);

end
