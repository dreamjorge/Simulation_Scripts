function [] = plotPropagatedRays(rayH11,rayH12,rayH21,RayH22)
% plot rays from hankels
  
  Nz = size(rayH11,2); %% number of points in z-direction

  
  switch nargin
    
    case 2 %%case of 2 hankels
      hold on
      for z_index = 1:Nz-1
        
        scatter3(rayH11(z_index).zCoordinate,rayH11(z_index).xCoordinate,rayH11(z_index).yCoordinate,20,'filled','MarkerFaceColor','r')
        scatter3(rayH12(z_index).zCoordinate,rayH12(z_index).xCoordinate,rayH12(z_index).yCoordinate,20,'filled','MarkerFaceColor','m')
        
      end
      hold off
      
    case 4 %% case of 4 hankels
      hold on
      for z_index = 1:Nz-1
        
        scatter3(rayH11(z_index).zCoordinate,rayH11(z_index).xCoordinate,rayH11(z_index).yCoordinate,20,'filled','MarkerFaceColor','r')
        scatter3(rayH12(z_index).zCoordinate,rayH12(z_index).xCoordinate,rayH12(z_index).yCoordinate,20,'filled','MarkerFaceColor','m')
        scatter3(rayH12(z_index).zCoordinate,rayH21(z_index).xCoordinate,rayH21(z_index).yCoordinate,20,'filled','MarkerFaceColor','y')
        scatter3(rayH12(z_index).zCoordinate,rayH22(z_index).xCoordinate,rayH22(z_index).yCoordinate,20,'filled','MarkerFaceColor','c')
        
        
      end
      hold off   

  end
end