function [] = plotRaysSquare(rayObject,color,scale)
%% Function for plot rays in a squre, usefull for Hermite Rays Propagatation 
  hold on

  plot([rayObject.xCoordinate/scale,rayObject.xCoordinate(1)/scale]...
      ,[rayObject.yCoordinate/scale,rayObject.yCoordinate(1)/scale]...
      ,'-.','LineWidth',3,'Color',color ...
      )

  pause(0.05)
  hold off
end