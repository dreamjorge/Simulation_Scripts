function h=plotOpticalField(xAxis,yAxis,OpticalField,ColorMap,units)
%% funtion that plots Optical Field xAxis,yAxis,OpticalField,ColorMap,units as inputs.
% Inputs:
%  - xAxis        as vector
%  - xAxis        as vector
%  - OpticalField as matrix
%  - ColorMap     as matrix
%  - units        as string 

  h = imagesc(xAxis, yAxis,OpticalField);
  colormap(ColorMap)
  axis square
  %shading flat
  labelxaxis=['$x$'];
  xlabel(labelxaxis,'Interpreter','latex','FontSize',18)
  labelyaxis=['$y$'];
  ylabel(labelyaxis,'Interpreter','latex','FontSize',18)
  set(gca,'YDir','normal')
  set(get(gca,'YLabel'),'Rotation',0)
end 