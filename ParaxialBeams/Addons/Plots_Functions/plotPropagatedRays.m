function [] = plotPropagatedRays(rayH11,rayH12,rayH21,rayH22,scale1,scale2,z_index)
% plot rays from hankels
  
  Nz = size(rayH11,2); %% number of points in z-direction
  Nrays = size(rayH11(1).zCoordinate,2);                      

  
  
  switch nargin
    
    case 2 %%case of 2 hankels
      hold on
      for z_index = 1:Nz-1
        
        scatter3(rayH11(z_index).zCoordinate,rayH11(z_index).xCoordinate,rayH11(z_index).yCoordinate,20,'filled','MarkerFaceColor','r')
        scatter3(rayH12(z_index).zCoordinate,rayH12(z_index).xCoordinate,rayH12(z_index).yCoordinate,20,'filled','MarkerFaceColor','m')
        
      end
      hold off
      
    case 6 %% case of 4 hankels
      
      
      hold on
      for z_index = 1:Nz
        
%         rayH11.(z_index).zCoordinate(Nrays+1) = rayH11.(z_index).zCoordinate(1);
%         rayH21.(z_index).zCoordinate(Nrays+1) = rayH21.(z_index).zCoordinate(1);
%         rayH12.(z_index).zCoordinate(Nrays+1) = rayH12.(z_index).zCoordinate(1);
%         rayH22.(z_index).zCoordinate(Nrays+1) = rayH22.(z_index).zCoordinate(1);
%         
        
        s1 = plot3([rayH11(z_index).zCoordinate/scale2,rayH11(z_index).zCoordinate(1)/scale2]...
                  ,[rayH11(z_index).xCoordinate/scale1,rayH11(z_index).xCoordinate(1)/scale1]...
                  ,[rayH11(z_index).yCoordinate/scale1,rayH11(z_index).yCoordinate(1)/scale1]...
                  ,'r');
        
        s2 = plot3([rayH12(z_index).zCoordinate/scale2,rayH12(z_index).zCoordinate(1)/scale2]...
                  ,[rayH12(z_index).xCoordinate/scale1,rayH12(z_index).xCoordinate(1)/scale1]...
                  ,[rayH12(z_index).yCoordinate/scale1,rayH12(z_index).yCoordinate(1)/scale1]...        
                  ,'m');
        
        s3 = plot3([rayH21(z_index).zCoordinate/scale2,rayH21(z_index).zCoordinate(1)/scale2]...
                  ,[rayH21(z_index).xCoordinate/scale1,rayH21(z_index).xCoordinate(1)/scale1]...
                  ,[rayH21(z_index).yCoordinate/scale1,rayH21(z_index).yCoordinate(1)/scale1]...        
                  ,'y');

        s4 = plot3([rayH22(z_index).zCoordinate/scale2,rayH22(z_index).zCoordinate(1)/scale2]...
                  ,[rayH22(z_index).xCoordinate/scale1,rayH22(z_index).xCoordinate(1)/scale1]...
                  ,[rayH22(z_index).yCoordinate/scale1,rayH22(z_index).yCoordinate(1)/scale1]...        
                  ,'c');                
        
      end
      hold off   

  end
end