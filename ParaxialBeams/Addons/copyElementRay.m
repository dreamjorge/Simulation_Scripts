function [rayOuput]= copyElementRay(rayInput,...
                                    index_point)
%% function for copy ray object
  
 
  rayOuput.xCoordinate    (index_point) = rayInput.xCoordinate    (index_point);
  rayOuput.yCoordinate    (index_point) = rayInput.yCoordinate    (index_point);
  rayOuput.zCoordinate    (index_point) = rayInput.zCoordinate    (index_point);
  rayOuput.rCoordinate    (index_point) = rayInput.rCoordinate    (index_point);
  rayOuput.thetaCoordinate(index_point) = rayInput.thetaCoordinate(index_point);
  
  if (index_point>1)
    rayOuput.zxSlope        (index_point) = rayInput.zxSlope        (index_point);
    rayOuput.zySlope        (index_point) = rayInput.zxSlope        (index_point);
    rayOuput.xySlope        (index_point) = rayInput.zxSlope        (index_point);
  else
    rayOuput.zxSlope        (index_point) = Inf;
    rayOuput.zySlope        (index_point) = Inf;
    rayOuput.xySlope        (index_point) = Inf;
  end
end