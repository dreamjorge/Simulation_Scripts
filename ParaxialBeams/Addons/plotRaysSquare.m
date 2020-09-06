function [] = plotRaysSquare(rayObject,color,scale)

hold on


plot([rayObject.xCoordinate/scale,rayObject.xCoordinate(1)/scale]...
    ,[rayObject.yCoordinate/scale,rayObject.yCoordinate(1)/scale]...
    ,'LineWidth',3,'Color',color ...
    )

%scatter(rayObject.xCoordinate/scale,rayObject.yCoordinate/scale,150,color,'filled','o')
pause(0.05)
hold off

end