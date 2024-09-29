function h = plotOpticalField(xAxis,yAxis,OpticalField,ColorMap,xLabelInput,yLabelInput)
%% funtion that plots Optical Field
%  it helps of imagesc, includes correction on axis, and input of colormap
%  and labels on x and y.
%  It uses xAxis,yAxis,OpticalField,ColorMap,xlabel,ylabel as inputs.
%  It returns imagesc object.
%
% Inputs:
%  - xAxis        as vector
%  - xAxis        as vector
%  - OpticalField as matrix
%  - ColorMap     as matrix
%  - xlabel       as string, we include latex interpeter
%  - ylabel       as string, we include latex interpeter

  h = imagesc(xAxis, yAxis,OpticalField);
  colormap(ColorMap)
  xlabel(xLabelInput,'Interpreter','latex','FontSize',18)
  ylabel(yLabelInput,'Interpreter','latex','FontSize',18)
  set(gca,'YDir','normal')
%   set(get(gca,'YLabel'),'Rotation',0) if we want rote y axis
end 