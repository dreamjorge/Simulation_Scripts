function    [Ray]=assignCoordinates2Ray(xi,yi,zi,Ray,point_index)
% Obtain Optical Ray with (xi,yi,zi) coordinates 

  %rayTotal = numel(Ray.xCoordinate);

  Ray.xCoordinate(point_index) = xi;
  Ray.yCoordinate(point_index) = yi;
  Ray.zCoordinate(point_index) = zi;
  
  [Ray.thetaCoordinate(point_index) ,Ray.rCoordinate(point_index)] = cart2pol(xi,yi); 
  
  %initial conditions of slopes
  Ray.xySlope(point_index) = Inf;
  Ray.zxSlope(point_index) = Inf;
  Ray.zySlope(point_index) = Inf;
  
  
  %initial conditions of slopes
  Ray.rthSlope(point_index) = Inf;
  Ray.zrSlope(point_index)  = Inf;
  Ray.zthSlope(point_index) = Inf;
  
end
  