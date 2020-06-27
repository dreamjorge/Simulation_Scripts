function [] = plotRaysPropagated(rayH1,rayH2,Nz)

  for z_index = 1:Nz-1
hold on

    scatter3(rayH1(z_index).zCoordinate,rayH1(z_index).xCoordinate,rayH1(z_index).yCoordinate,20,'filled','MarkerFaceColor',[0 .75 .75])
    scatter3(rayH2(z_index).zCoordinate,rayH2(z_index).xCoordinate,rayH2(z_index).yCoordinate,20,'filled','MarkerFaceColor',[1 0 0])

  end
  hold off

  xlabel('z[microns]')
  ylabel('x[microns]')
  zlabel('y[microns]')

end