function    [Ray]=assignCoordinates2CylindricalRay(xi,yi,zi,Ray,point_index,hankelType)

  Ray.xCoordinate(point_index) = xi;
  Ray.yCoordinate(point_index) = yi;
  Ray.zCoordinate(point_index) = zi;
  
  %initial conditions of slopes
  Ray.rthSlope(point_index) = Inf;
  Ray.zrSlope(point_index)  = Inf;
  Ray.zthSlope(point_index) = Inf;
  
  % hankel type of rays
  Ray.hankelType(point_index) = hankelType;
  
end
  