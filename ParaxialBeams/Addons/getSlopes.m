function [ray] = getSlopes(ray,x,y,z,...
                               dx,dy,dx,...
                               xi,yi,zi,...
                               InitialWaist,Wavelength,l,p,nh)

%y = cte, z = cte
HLx = HankelLaguerre(x,yi,zi,InitialWaist,Wavelength,l,p,nh);

%x = cte, z = cte
HLy = HankelLaguerre(xi,y,zi,InitialWaist,Wavelength,l,p,nh);

%x = cte, y = cte
HLz = HankelLaguerre(xi,yi,z,InitialWaist,Wavelength,l,p,nh);

fx   = unwrap(angle(HLx));
fy   = unwrap(angle(HLy));
fz   = unwrap(angle(HLz));

[ray.zxSlope,ray(jj).zySlope,ray(jj).xySlope] = gradientxyz(fx,fy,fz,k,dx,dx,dz,xi,yi,xi);

end
