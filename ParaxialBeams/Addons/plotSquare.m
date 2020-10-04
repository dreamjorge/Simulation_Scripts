function [] = plotSquare(x,y,color)

hold on

plot([x,x(1)]...
    ,[y,y(1)]...
    ,'LineWidth',3,'Color',color ...
    )
hold off

end