function [rayObject] = getLaguerreSlopes(rayObject,x,y,z,...
                                         dx,dy,dz,...
                                         xi,yi,zi,...
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
                                                 
  %y = cte, z = cte
  HLx = HankelLaguerre(x,yi,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
  
  %x = cte, z = cte
  HLy = HankelLaguerre(xi,y,zi,InitialWaist,Wavelength,nu,mu,numberHankel);

  %x = cte, y = cte
  HLz = HankelLaguerre(xi,yi,z,InitialWaist,Wavelength,nu,mu,numberHankel);

  fx  = unwrap(angle(HLx.OpticalField));
  fy  = unwrap(angle(HLy.OpticalField));
  fz  = unwrap(angle(HLz.OpticalField));

  %Estimate k with wavelength
  k   = (2*pi)/Wavelength;
  
  %Calculating gradient for slopes of vector normal
  [rayObject.zxSlope, rayObject.zySlope, rayObject.xySlope] = ...
   gradientxyz(fx,fy,fz,k,dx,dy,dz,xi,yi,zi);

end
