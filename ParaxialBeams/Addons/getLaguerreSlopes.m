function [Slopes] = getLaguerreSlopes(rayObject,...
                                      r,th,z,...
                                      dr,...
                                      LaguerreParametersObject,...
                                      numberHankel)
                                       
% This function estimates components (x,y,z) of a ray object and its slopes
% for next iteration.
%
% receives  rayObject:                Ray object where includes (x,y,z) and (r,th,z ) 
%                                     coordinates.
%           r:                        Distance of point to center of coordinates vector 
%                                     of coordinates.
%           th:                       Angle from straight line from point to center of 
%                                     coordinates to horizontal axis.
%           z:                        Distance of z coordinate to center of coordinates.   
%           dr:                       diferential vector (dx,dy,dz)
%           LaguerreParametersObject: Paramerters of Laguerre calculated for zi.
%           numberHankel              Type Hankel (1 of 2).
% gives     rayObject properties calculated.
%
% Example : [rayObject] = getLaguerreSlopes(rayObject,
%                                           r,th,z,...
%                                           dr,...
%                                           LaguerreParametersObject,...
%                                           numberHankel)
 

  %taking coordinates of point to calculate slope
  ri  = rayObject.rCoordinate;
  thi = rayObject.thetaCoordinate;
  %zi  = rayObject.zCoordinate;
  
  %transform coordinates to cartesians
  [xi,yi] = pol2cart(thi,ri);
  [x,y]   = pol2cart(th,r);
  
  %Parameter calculated to zi was calculated outside, and this is an input for this algo.
  LPzi    = LaguerreParametersObject;
  
  %extrating constants of next parameters
  InitialWaist = LPzi.InitialWaist;
  Wavelength   = LPzi.Wavelength;
  l            = LPzi.l;
  p            = LPzi.p;
  k            = LPzi.k;
  
  %Calculating parameters to z vector.
  LPz     = LaguerreParameters(z,InitialWaist,Wavelength,l,p);
  
  HLx     = HankelLaguerre(sqrt(x.^2 + yi.^2),atan2(yi,x) ,LPzi,numberHankel);
  HLy     = HankelLaguerre(sqrt(xi.^2+ y.^2) ,atan2(y,xi) ,LPzi,numberHankel);
  
  HLz     = HankelLaguerre(ri                ,thi         ,LPz ,numberHankel);
  
  fx      = unwrap(angle(HLx.OpticalFieldLaguerre));
  fy      = unwrap(angle(HLy.OpticalFieldLaguerre));
  fz      = unwrap(angle(HLz.OpticalFieldLaguerre));

 
  [Slopes] = gradientxyz(fx,fy,fz,k,dr,rayObject);
  
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
