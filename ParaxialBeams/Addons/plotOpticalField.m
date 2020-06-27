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
  labelxaxis=['$x \left[ ',units,' \right]$'];
  xlabel(labelxaxis,'Interpreter','latex')
  labelyaxis=['$y \left[ ',units,' \right]$'];
  ylabel(labelyaxis,'Interpreter','latex') 
  set(gca,'YDir','normal')
end 