function [h] = plotPropagatedRays3D(rayObject,scaleX,scaleY,scaleZ,color)
% plot all rays in 3D  

  Nrays = size(rayObject(1).zCoordinate,2);
  Nz    = size(rayObject,2); %% number of points in z-direction
  ray   = struct();

  for z_index = 1:Nz

    for ray_index = 1:Nrays
          
      ray(ray_index).z(z_index) = rayObject(z_index).zCoordinate(ray_index);
      ray(ray_index).x(z_index) = rayObject(z_index).xCoordinate(ray_index);
      ray(ray_index).y(z_index) = rayObject(z_index).yCoordinate(ray_index);
          
    end
 
  end  
  
  hold on
      
  for ray_index = 1:Nrays
        
    h = plot3(scaleZ*ray(ray_index).z,...
              scaleX*ray(ray_index).x,...
              scaleY*ray(ray_index).y,...
              'Color',color,'LineWidth',2);
      
  end
      
  hold off
  
end