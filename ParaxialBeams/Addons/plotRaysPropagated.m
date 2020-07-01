function [] = plotPropagatedRays(rayH11,rayH12,rayH21,RayH22)
% plot rays

  switch nargin
    
    case 3 
      hold on
      for z_index = 1:Nz-1
        
        scatter3(rayH11(z_index).zCoordinate,rayH11(z_index).xCoordinate,rayH11(z_index).yCoordinate,20,'filled','MarkerFaceColor',[0 .75 .75])
        scatter3(rayH12(z_index).zCoordinate,rayH12(z_index).xCoordinate,rayH12(z_index).yCoordinate,20,'filled','MarkerFaceColor',[1 0 0])
        
      end
      hold off

end