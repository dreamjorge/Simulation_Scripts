function h = plotCircleIn3D(x,y,z,r,color,lineWidth)
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
  zunit = zeros(1,numel(th));
  h = plot3(xunit, yunit,z+zunit,'Color',color,'LineWidth',lineWidth);
  hold off
end