function [] = plotRays(rayObject,color)

hold on

scatter(rayObject.xCoordinate,rayObject.yCoordinate,10,color,'filled','o')
pause(0.05)
hold off

end