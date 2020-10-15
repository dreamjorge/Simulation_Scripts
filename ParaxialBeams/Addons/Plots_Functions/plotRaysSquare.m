function [h] = plotRaysSquare(rayObject,color,scaleX,scaleY)
%% Function for plot rays in a squre, usefull for Hermite Rays Propagatation 
  hold on

  h = plot([rayObject.xCoordinate*scaleX,rayObject.xCoordinate(1)*scaleY]...
          ,[rayObject.yCoordinate*scaleX,rayObject.yCoordinate(1)*scaleY]...
          ,'-.','LineWidth',3,'Color',color ...
          );

  pause(0.05)
  hold off
end