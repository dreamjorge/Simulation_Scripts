function [rayObject] = propagateRay(rayObject,actualStep,dz)

  previouStep = actualStep-1;
  
  rayObject.thetaCoordinate(actualStep) = rayObject.thetaCoordinate(previouStep) + (1./rayObject.zthSlope)*dz;
  rayObject.rCoordinate(actualStep)     = rayObject.rCoordinate(previouStep)     + (1./rayObject.zrSlope) *dz;
  

end