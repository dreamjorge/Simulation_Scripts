function [s] = plotPropagatedRaysSquareInZfor3D(rayObject,scaleX,scaleY,scaleZ,z_index,color)
%% plot rays in square from hankels in 3D

hold on

  x = [rayObject(z_index).xCoordinate,rayObject(z_index).xCoordinate(1)];
  y = [rayObject(z_index).yCoordinate,rayObject(z_index).yCoordinate(1)];
  z = [rayObject(z_index).zCoordinate,rayObject(z_index).zCoordinate(1)];

  s = plot3(scaleZ*z,scaleX*x,scaleY*y,'Color',color,'LineWidth',2);

hold off

end