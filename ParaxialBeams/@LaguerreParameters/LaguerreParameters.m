classdef LaguerreParameters < GaussianParameters
% This class gives Laguerre parameters of Laguerre Gaussian Beam
% recives PropagationDistance,InitialWaist,Wavelength,l,p as input and
% gives a object with next properties:
%
% -l
% -p
% -LaguerreWaist
% -PhiPhase
% -PropagationDistance
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
% Example: Parameters = GaussianParameters(PropagationDistance,InitialWaist,Wavelength)
  
  
  properties
    l
    p
  end
  
  properties (Dependent)
    LaguerreWaist
    PhiPhase
  end
  
  methods(Static)
    waistL = waistFunction(PropagationDistance,InitialWaist,RayleighDistance,nu,mu);
  end
  
  
  methods
    
    function PhiPhase = get.PhiPhase(obj) 
      PhiPhase = (abs(obj.l)+2*obj.p).*obj.GouyPhase;
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = LaguerreParameters.waistFunction(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end

    %% Constructor
    function Parameters = LaguerreParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu)

     Parameters@GaussianParameters(PropagationDistance,InitialWaist,Wavelength);

     if nargin == 5 
       Parameters.l             = nu;
       Parameters.p             = mu;
     else
        error('You need introduce PropagationDistance, InitialWaist and Wavelength, l, p inputs')
     end

    end
    
   end
end