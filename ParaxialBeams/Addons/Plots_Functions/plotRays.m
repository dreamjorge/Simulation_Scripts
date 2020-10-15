function [h] = plotRays(rayObject,color,scaleX,scaleY)
%% Function for plot rays in x,y plane
  hold on

  h = scatter(rayObject.xCoordinate*scaleX,rayObject.yCoordinate*scaleY,10,color,'filled','o');

  hold off

end