function [] = plotPropagatedRays2inZ(rayH11,rayH12,rayH21,rayH22,z_index)
  Nz    = size(rayH11,2); %% number of points in z-direction
  Nrays = size(rayH11(1).zCoordinate,2);                      

  
  ray11 = struct();
  ray12 = struct();
  ray21 = struct();
  ray22 = struct();
  for z_index = 1:Nz

        for ray_index = 1:Nrays
          
          ray11(ray_index).z(z_index) = rayH11(z_index).zCoordinate(ray_index);
          ray11(ray_index).x(z_index) = rayH11(z_index).xCoordinate(ray_index);
          ray11(ray_index).y(z_index) = rayH11(z_index).yCoordinate(ray_index);
          
          ray12(ray_index).z(z_index) = rayH12(z_index).zCoordinate(ray_index);
          ray12(ray_index).x(z_index) = rayH12(z_index).xCoordinate(ray_index);
          ray12(ray_index).y(z_index) = rayH12(z_index).yCoordinate(ray_index);
          
          ray21(ray_index).z(z_index) = rayH21(z_index).zCoordinate(ray_index);
          ray21(ray_index).x(z_index) = rayH21(z_index).xCoordinate(ray_index);
          ray21(ray_index).y(z_index) = rayH21(z_index).yCoordinate(ray_index);
          
          ray22(ray_index).z(z_index) = rayH22(z_index).zCoordinate(ray_index);
          ray22(ray_index).x(z_index) = rayH22(z_index).xCoordinate(ray_index);
          ray22(ray_index).y(z_index) = rayH22(z_index).yCoordinate(ray_index);
          
        end
 
  end    
      hold on
      
      for ray_index = 1:Nrays
        
       %  plot3(ray11(ray_index).z,ray11(ray_index).x,ray11(ray_index).y,'r','LineWidth',2)
         plot3(ray12(ray_index).z,ray12(ray_index).x,ray12(ray_index).y,'m','LineWidth',2)
       %  plot3(ray21(ray_index).z,ray21(ray_index).x,ray21(ray_index).y,'y','LineWidth',2)
        % plot3(ray22(ray_index).z,ray22(ray_index).x,ray22(ray_index).y,'c','LineWidth',2)
        
      end
      
      hold off
end