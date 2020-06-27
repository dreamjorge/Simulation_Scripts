function [Slopes] = getLaguerreSlopes(rayObject,...
                                      x,y,... transversal coordinates
                                      dr,...
                                      LaguerreParametersZiObject,...
                                      LaguerreParametersZObject,...
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
 


  
  fx      = unwrap(angle(HLx.OpticalFieldLaguerre));
  fy      = unwrap(angle(HLy.OpticalFieldLaguerre));
  fz      = unwrap(angle(HLz.OpticalFieldLaguerre));

 
  [rayObject] = gradientxyz(fx,fy,fz,LaguerreParametersZiObject.k,dr,rayObject);
  
end
