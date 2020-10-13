function h = plotCircle(x,y,r,color,lineWidth)
%% function recives 
% x         - center position of circle to draw
% y         - center position of circle to draw
% r         - radius of circle to draw
% color     - color to draw
% lineWidth - width of lone to draw circle
  hold on
  th = 0:pi/50:2*pi;
  xunit = r * cos(th) + x;
  yunit = r * sin(th) + y;
  h = plot(xunit, yunit,'Color',color,'LineWidth',lineWidth);
  hold off
end