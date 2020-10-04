function [] = plotPropagatedRaysInZfor3D(rayH11,rayH12,rayH21,rayH22,scale1,scale2,z_index)
% plot rays from hankels
  
  Nz = size(rayH11,2); %% number of points in z-direction
  Nrays = size(rayH11(1).zCoordinate,2);                      

  
  
  switch nargin
    
    case 2 %%case of 2 hankels
      hold on
      for z_index = 1:Nz-1
        
        scatter3(rayH11(z_index).xCoordinate,rayH11(z_index).yCoordinate,20,'filled','MarkerFaceColor','r')
        scatter3(rayH12(z_index).zCoordinate,rayH12(z_index).xCoordinate,rayH12(z_index).yCoordinate,20,'filled','MarkerFaceColor','m')
        
      end
      hold off
      
    case 7 %% case of 4 hankels
      
      
      hold on
%       for z_index = 1:Nz
        
%         rayH11.(z_index).zCoordinate(Nrays+1) = rayH11.(z_index).zCoordinate(1);
%         rayH21.(z_index).zCoordinate(Nrays+1) = rayH21.(z_index).zCoordinate(1);
%         rayH12.(z_index).zCoordinate(Nrays+1) = rayH12.(z_index).zCoordinate(1);
%         rayH22.(z_index).zCoordinate(Nrays+1) = rayH22.(z_index).zCoordinate(1);
%         
        x = [rayH12(z_index).xCoordinate/scale1,rayH12(z_index).xCoordinate(1)/scale1];
        y = [rayH12(z_index).yCoordinate/scale1,rayH12(z_index).yCoordinate(1)/scale1];
        z = [rayH12(z_index).zCoordinate/scale1,rayH12(z_index).zCoordinate(1)/scale1];

        s1 = plot3(z,y,x,'m','LineWidth',2);
        
%         s2 = plot([rayH12(z_index).xCoordinate/scale1,rayH12(z_index).xCoordinate(1)/scale1]...
%                   ,[rayH12(z_index).yCoordinate/scale1,rayH12(z_index).yCoordinate(1)/scale1]...        
%                   ,'m','LineWidth',2);
%         
%         s3 = plot([rayH21(z_index).xCoordinate/scale1,rayH21(z_index).xCoordinate(1)/scale1]...
%                   ,[rayH21(z_index).yCoordinate/scale1,rayH21(z_index).yCoordinate(1)/scale1]...        
%                   ,'y','LineWidth',2);
% 
%         s4 = plot([rayH22(z_index).xCoordinate/scale1,rayH22(z_index).xCoordinate(1)/scale1]...
%                   ,[rayH22(z_index).yCoordinate/scale1,rayH22(z_index).yCoordinate(1)/scale1]...        
%                   ,'c','LineWidth',2);                
        
%       end
      hold off   

  end
end