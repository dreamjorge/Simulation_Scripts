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
  label=['$x \left[ ',convertCharsToStrings(units),' \right]$'];
  label=strjoin(label);
  xlabel(label,'Interpreter','latex')
  label=['$y \left[ ',convertCharsToStrings(units),' \right]$'];
  label=strjoin(label);
  ylabel(label,'Interpreter','latex') 
  set(gca,'YDir','normal')
end 