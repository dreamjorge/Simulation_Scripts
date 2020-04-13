function [rayObject] = getLaguerreSlopes(rayObject,...
                                         r,th,z,...
                                         dr,...
                                         InitialWaist,Wavelength,nu,mu,numberHankel)
                                       
% This function estimates components (x,y,z) of a ray object and its slopes
% for next iteration.
%
% receives  rayObject,x,y,z,dx,dy,dz,xi,yi,zi,InitialWaist,Wavelength,nu,mu,
%           and numberHankel as inputs.
% gives     rayObject properties calculated.
%
% Example : [rayObject] = getLaguerreSlopes(rayObject,x,y,z,...
%                                           dx,dy,dz,...
%                                           xi,yi,zi,...
%                                           InitialWaist,Wavelength,nu,mu,numberHankel)
       

  ri   = rayObject.rCoordinate;
  thi  = rayObject.thetaCoordinate;
  zi   = rayObject.zCoordinate;
  
  %transform coordinates to cartsians
  [xi,yi] = pol2cart(thi,ri);
  [x,y]   = pol2cart(th,r);
  
  %vriable
  
  HLx     = HankelLaguerre(sqrt(x.^2 + yi.^2),atan2(yi,x) ,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
  HLy     = HankelLaguerre(sqrt(xi.^2+ y.^2) ,atan2(y,xi) ,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
  HLz     = HankelLaguerre(sqrt(xi.^2+ yi.^2),atan2(yi,xi),z ,InitialWaist,Wavelength,nu,mu,numberHankel);
  
  fx      = unwrap(angle(HLx.OpticalField));
  fy      = unwrap(angle(HLy.OpticalField));
  fz      = unwrap(angle(HLz.OpticalField));

  k   = (2*pi)/Wavelength;
  
  [rayObject] = gradientxyz(fx,fy,fz,k,dr,rayObject);
  
%   %th = cte, z = cte
%   HLr  = HankelLaguerre(r,thi,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
%   
%   %r = cte, z = cte
%   HLth = HankelLaguerre(ri,th,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
% 
%   %r = cte, th = cte
%   HLz  = HankelLaguerre(ri,thi,z,InitialWaist,Wavelength,nu,mu,numberHankel);
% 
%   fr   = unwrap(angle(HLr.OpticalField));
%   fth  = unwrap(angle(HLth.OpticalField));
%   fz   = unwrap(angle(HLz.OpticalField));
% 
%   %Estimate k with wavelength
%   k   = (2*pi)/Wavelength;
  
%   %Calculating gradient for slopes of vector normal
%   [rayObject] = ...
%    gradientCylindrical(fr,fth,fz,k,dr,rayObject);

end
