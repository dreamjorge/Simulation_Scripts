function [] = plotRays(rayObject,color)

rayTotal = numel(rayObject.xCoordinate);

hold on
for point_index=1:rayTotal
    plot(rayObject.xCoordinate(point_index)...
        ,rayObject.yCoordinate(point_index)...
        ,'.','MarkerSize',10,'LineWidth',2,'color',color)
end
hold off

end