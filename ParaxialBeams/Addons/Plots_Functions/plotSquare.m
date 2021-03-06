function [h] = plotSquare(x,y,color)
%% Function for plot square given vertices in x,y vectors

  hold on

  h = plot([x,x(1)],[y,y(1)],'LineWidth',3,'Color',color);
  
  hold off

end