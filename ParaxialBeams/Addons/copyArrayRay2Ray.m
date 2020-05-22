function [rayOutput]= copyArrayRay2Ray(rayInput,...
                                       rayOutput,...
                                       index_point)
                            
                            
  rayOutput.xCoordinate     = rayInput.xCoordinate    (index_point);
  rayOutput.yCoordinate     = rayInput.yCoordinate    (index_point);
  rayOutput.zCoordinate     = rayInput.zCoordinate    (index_point);
  rayOutput.rCoordinate     = rayInput.rCoordinate    (index_point);
  rayOutput.thetaCoordinate = rayInput.thetaCoordinate(index_point);
  
  rayOutput.zxSlope  = rayInput.zxSlope(index_point);
  rayOutput.zySlope  = rayInput.zySlope(index_point);
  rayOutput.xySlope  = rayInput.xySlope(index_point);
  
  end