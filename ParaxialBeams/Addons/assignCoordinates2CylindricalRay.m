function    [Ray] = assignCoordinates2CylindricalRay(xi,yi,zi,Ray,point_index,hankelType)

  [thi,ri]  = cart2pol(xi,yi);


  Ray.rCoordinate(point_index)       =   ri;
  Ray.thetaCoordinate(point_index)   =  thi;

  Ray.xCoordinate(point_index) = xi;
  Ray.yCoordinate(point_index) = yi;
  Ray.zCoordinate(point_index) = zi;
  
  %initial conditions of slopes
  Ray.zrSlope(point_index)  = Inf;
  Ray.zthSlope(point_index) = Inf;
  Ray.rthSlope(point_index) = Inf;
  
  % hankel type of rays
  Ray.hankelType(point_index) = hankelType;
  
end
