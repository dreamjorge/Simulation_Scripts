function [ray] = getHermiteSlopes(ray,r,dr,ri,HermiteParameters,Hankel)

  nr_Hankelx = Hankel(1,1);
  nr_Hankely = Hankel(1,2);

  %y = cte, z = cte
  HLx = HankelHermite(r.x,ri.y,ri.z,HermiteParameters.InitialWaist,HermiteParameters.Wavelength,...
                      HermiteParameters.n,HermiteParameters.m,nr_Hankelx,nr_Hankely);
  
  %x = cte, z = cte
  HLy = HankelHermite(ri.x,r.y,ri.z,HermiteParameters.InitialWaist,HermiteParameters.Wavelength,...
                      HermiteParameters.n,HermiteParameters.m,nr_Hankelx,nr_Hankely);
  
  %x = cte, y = cte
  HLz = HankelHermite(ri.x,ri.y,r.z,HermiteParameters.InitialWaist,HermiteParameters.Wavelength,...
                      HermiteParameters.n,HermiteParameters.m,nr_Hankelx,nr_Hankely);
  
  fx  = unwrap(angle(HLx.OpticalField));
  fy  = unwrap(angle(HLy.OpticalField));
  fz  = unwrap(angle(HLz.OpticalField));

  %k   = (2*pi)/HermiteParameters.InitialWaist;
  
  [ray.zxSlope,ray.zySlope,ray.xySlope] = gradientxyz(fx,fy,fz,HermiteParameters.k,dr.x,dr.y,dr.z,ri.x,ri.y,ri.z);

end