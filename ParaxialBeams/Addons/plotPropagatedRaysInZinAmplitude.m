function [] = plotPropagatedRaysInZinAmplitude(rayH11,rayH12,rayH21,rayH22,scale1,scale2,z_index)

  
  Nz = size(rayH11,2); %% number of points in z-direction
  Nrays = size(rayH11(1).zCoordinate,2);                      

  
        hold on
 
        s1 = plot([rayH11(z_index).xCoordinate/scale1,rayH11(z_index).xCoordinate(1)/scale1]...
                  ,[rayH11(z_index).yCoordinate/scale1,rayH11(z_index).yCoordinate(1)/scale1]...
                  ,'r','LineWidth',2);
        
        s2 = plot([rayH12(z_index).xCoordinate/scale1,rayH12(z_index).xCoordinate(1)/scale1]...
                  ,[rayH12(z_index).yCoordinate/scale1,rayH12(z_index).yCoordinate(1)/scale1]...        
                  ,'m','LineWidth',2);
        
        s3 = plot([rayH21(z_index).xCoordinate/scale1,rayH21(z_index).xCoordinate(1)/scale1]...
                  ,[rayH21(z_index).yCoordinate/scale1,rayH21(z_index).yCoordinate(1)/scale1]...        
                  ,'y','LineWidth',2);

        s4 = plot([rayH22(z_index).xCoordinate/scale1,rayH22(z_index).xCoordinate(1)/scale1]...
                  ,[rayH22(z_index).yCoordinate/scale1,rayH22(z_index).yCoordinate(1)/scale1]...        
                  ,'c','LineWidth',2);                
        
%       end
      hold off   

  
  
  
  
  
end