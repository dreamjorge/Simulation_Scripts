function [h] = plotWaistHermite2D(HermiteParametersZi,scaleX,scaleY,color)
%% Function for plot waist of Hermite Gauss Beam in 2D 
% at z - distance given in parameters.
% Recives Hermite Parameters and returns plot object.
  
  % Parameters needed from Hermite Parameters
  rd     = HermiteParametersZi.RayleighDistance;
  wo     = HermiteParametersZi.InitialWaist;
  n      = HermiteParametersZi.n;
  m      = HermiteParametersZi.m;
  zi     = HermiteParametersZi.zCoordinate;
  % Waist 1-dimentional in each axes
  wx     = HermiteParametersZi.getWaistOneDirection(zi,wo,rd,n);
  wy     = HermiteParametersZi.getWaistOneDirection(zi,wo,rd,m);
  % Vertices of square
  xwaist = [-wx/2,  wx/2, wx/2, -wx/2];
  ywaist = [-wy/2, -wy/2, wy/2,  wy/2];
  % Plot Waist of Hermite Gauss Beam
  h      = plotSquare(xwaist*scaleX,ywaist*scaleY,color);
  
end