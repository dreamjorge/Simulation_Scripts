function [h] = plotPropagatedRayZXPlane(RayObject,ray_index,scale1,scale2,color,linewidth)
%% Function for plot propagated ray given by ray_index.
% inputs are:
% -RayObject    Object that includes all rays information propagated.
% -ray_index    Which ray we want in RayObject.
% -scale1       Scale for x-axis.
% -scale2       Scale for z-axis.
% -color        Color to plot.
% -linewidth    width of line to plot.
% 
% Example:
% plotPropagatedRayZYPlane(RayObject,ray_index,scale1,scale2,color).

  % number of points in z-direction
  Nz = size(RayObject,2);

  % Extract information of ray
  for z_index = 1:Nz
          
    rayTemp.z(z_index) = RayObject(z_index).zCoordinate(ray_index);
    rayTemp.x(z_index) = RayObject(z_index).xCoordinate(ray_index);
          
  end

  % Plot
  h = plot(rayTemp.z/scale2,rayTemp.x/scale1,...
          'color',color,'LineWidth',linewidth);

end