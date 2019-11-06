function [ray] = getHermiteSlopes(ray,r,dr,ri,HermiteParameters,Hankel)

nr_Hankelx = Hankel(1,1);

nr_Hankely = Hankel(1,2);

nargin

  %y = cte, z = cte
  HLx = HankelHermite(r.x,ri.y,ri.z,HermiteParameters.InitialWaist,HermiteParameters.Wavelength,...
                      HermiteParameters.n,HermiteParameters.m,nr_Hankelx,nr_Hankely);
  
  %x = cte, z = cte
  HLy = HankelHermite(xi,y,zi,InitialWaist,Wavelength,nu,mu,nr_Hankelx,nr_Hankely);

  %x = cte, y = cte
  HLz = HankelHermite(xi,yi,z,InitialWaist,Wavelength,nu,mu,nr_Hankelx,nr_Hankely);

  fx   = unwrap(angle(HLx.OpticalField));
  fy   = unwrap(angle(HLy.OpticalField));
  fz   = unwrap(angle(HLz.OpticalField));

  k = 2*pi/InitialWaist;
  
  [ray.zxSlope,ray.zySlope,ray.xySlope] = gradientxyz(fx,fy,fz,k,dr.dx,dy,dz,xi,yi,zi);

end
