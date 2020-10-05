classdef LaguerreParameters < GaussianParameters & handle & matlab.mixin.Copyable
% LaguerreParameters class gives Laguerre parameters of Laguerre Gaussian Beam
% recives zCoordinate,InitialWaist,Wavelength,l,p as input and
% gives a object with next properties:
%
% -l
% -p
% -LaguerreWaist
% -PhiPhase
% -zCoordinate
% -InitialWaist
% -Wavelength
% -RayleighDistance
% -k
% -Waist
% -Radius
% -GouyPhase
% -Amplitude
% -DivergenceAngle
% 

% Example: Parameters = LaguerreParameters(zCoordinate...
%                                         ,InitialWaist...
%                                         ,Wavelength)
  
  properties
    l % angular number
    p % radial number
  end
  
  properties (Dependent)
    LaguerreWaist  % Waist of Laguerre Gauss Beam
    PhiPhase       % Phase of Laguerre Gauss Beam
  end
  
  methods(Static)   
    waistL = getWaist(zCoordinate,InitialWaist,RayleighDistance,nu,mu); % Function for estimate wast of Laguerre Gauss Beam
  end
  
  
  methods
    
    function PhiPhase      = get.PhiPhase(obj) 
      % Function for estimate phase of Laguerre Gaussian Beam
      PhiPhase = (abs(obj.l)+2*(obj.p)-1).*obj.GouyPhase; 
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      % Function for estimate phase of Laguerre Gaussian Beam
      LaguerreWaist = LaguerreParameters.getWaist(obj.zCoordinate,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end

    
    
    
    function Parameters = LaguerreParameters(zCoordinate,InitialWaist,Wavelength,l,p)
    % Constructor of Laguerre Parameters Object 
     %Call Gaussian Parameters
     Parameters@GaussianParameters(zCoordinate,InitialWaist,Wavelength);

     if nargin == 5 
       % Add parameter of input to object
       Parameters.l             = l;
       Parameters.p             = p;
     else
        error('You need introduce zCoordinate, InitialWaist and Wavelength, l, p inputs')
     end

    end
    
   end
end