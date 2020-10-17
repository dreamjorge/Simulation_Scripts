function [h] = plotRaysAtZ(rayObject,scaleX,scaleY,lineWidth,color)
%% function for plot rays at z-distance in x,y plane
% recives next elements
% - rayObect    structure includes information of ray to plot.
% - lineWidth   width of line to plot in this case size of circle on
%               scatter.
% - color       color that we want plot.
% - scale       scale for rescale elements in ray.

  hold on
  
  h = scatter(scaleX*rayObject.xCoordinate,scaleY*rayObject.yCoordinate,lineWidth,color,'filled','o');

  hold off

end