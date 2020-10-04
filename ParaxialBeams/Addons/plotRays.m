function [] = plotRays(rayObject,color,scale)

hold on

scatter(rayObject.xCoordinate/scale,rayObject.yCoordinate/scale,10,color,'filled','o')
pause(0.05)
hold off

end