function [ray] = getHermiteSlopes(ray,x,y,z,...
                                  dx,dy,dz,...
                                  xi,yi,zi,...
                                  InitialWaist,Wavelength,nu,mu,nrHankelx,nrHankely)

  %y = cte, z = cte
  HLx = HankelHermite(x,yi,zi,InitialWaist,Wavelength,nu,mu,nrHankelx,nrHankely);
  
  %x = cte, z = cte
  HLy = HankelHermite(xi,y,zi,InitialWaist,Wavelength,nu,mu,nrHankelx,nrHankely);

  %x = cte, y = cte
  HLz = HankelHermite(xi,yi,z,InitialWaist,Wavelength,nu,mu,nrHankelx,nrHankely);

  fx   = unwrap(angle(HLx.OpticalField));
  fy   = unwrap(angle(HLy.OpticalField));
  fz   = unwrap(angle(HLz.OpticalField));

  k = 2*pi/InitialWaist;
  
  [ray.zxSlope,ray.zySlope,ray.xySlope] = gradientxyz(fx,fy,fz,k,dx,dy,dz,xi,yi,zi);

end
