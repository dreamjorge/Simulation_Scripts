function [rayObject] = getNewRayCoordinates(rayObject,dr)

  rayObject.thetaCoordinate(ii) = rayObject.thetaCoordinate(ii-1) + (1./rayObject.zthSlope)*dz;
  rayObject.rCoordinate(ii)     = rayObject.rCoordinate(ii-1)     + (1./rayObject.zrSlope)*dz;

end