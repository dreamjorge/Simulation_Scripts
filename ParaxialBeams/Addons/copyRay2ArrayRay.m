function [rayOutput]= copyRay2ArrayRay(rayInput,...
                                       rayOutput,...
                                       index_output_point)
                                             
  rayOutput.xCoordinate    (index_output_point) = rayInput.xCoordinate    ;
  rayOutput.yCoordinate    (index_output_point) = rayInput.yCoordinate    ;
  rayOutput.zCoordinate    (index_output_point) = rayInput.zCoordinate    ;
  rayOutput.rCoordinate    (index_output_point) = rayInput.rCoordinate    ;
  rayOutput.thetaCoordinate(index_output_point) = rayInput.thetaCoordinate;

  rayOutput.zxSlope(index_output_point) = rayInput.zxSlope;
  rayOutput.zySlope(index_output_point) = rayInput.zySlope;
  rayOutput.xySlope(index_output_point) = rayInput.xySlope;

end