function [rayObject] = getLaguerreSlopes(rayObject,r,th,z,...
                                         dr,dth,dz,...
                                         ri,thi,zi,...
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
                                                 
  %th = cte, z = cte
  HLr = HankelLaguerre(r,thi,zi,InitialWaist,Wavelength,nu,mu,numberHankel);
  
  %x = cte, z = cte
  HLth = HankelLaguerre(ri,th,zi,InitialWaist,Wavelength,nu,mu,numberHankel);

  %x = cte, y = cte
  HLz = HankelLaguerre(ri,thi,z,InitialWaist,Wavelength,nu,mu,numberHankel);

  fr   = unwrap(angle(HLr.OpticalField));
  fth  = unwrap(angle(HLth.OpticalField));
  fz   = unwrap(angle(HLz.OpticalField));

  %Estimate k with wavelength
  k   = (2*pi)/Wavelength;
  
  %Calculating gradient for slopes of vector normal
  [rayObject.zrSlope, rayObject.zthSlope, rayObject.rthSlope] = ...
   gradientCylindrical(fr,fth,fz,k,dr,dth,dz,ri,thi,zi);

end
