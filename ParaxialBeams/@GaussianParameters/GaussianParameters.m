classdef GaussianParameters <  handle
% This class gives Gaussian parameters of Gaussian Beam
% recives PropagationDistance,InitialWaist,Wavelength as input and
% gives a object with next properties:
%
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
    PropagationDistance
    InitialWaist
    Wavelength
  end
  
  properties (Dependent)
    RayleighDistance
    k
    Waist
    Radius
    GouyPhase
    Amplitude
    DivergenceAngle
  end

  
  methods(Static)
%     
        Waist  = waistFunction(PropagationDistance,...
                               InitialWaist,...
                               RayleighDistance);
                                    
        Phase  = phaseFunction(PropagationDistance,...
                               RayleighDistance);
                                         
        Radius = radiusFunction(PropagationDistance,...
                                RayleighDistance);

  end
  
  
  methods
    
  %% Constructor
    function Parameters = GaussianParameters(PropagationDistance,InitialWaist,Wavelength)
      if nargin == 3 
        Parameters.PropagationDistance = PropagationDistance;
        Parameters.InitialWaist        = InitialWaist;
        Parameters.Wavelength          = Wavelength;
      else
        error('You need introduce PropagationDistance, InitialWaist and Wavelength inputs')
      end
    end

    

    %% Error for dependent properties
    
    function obj = set.RayleighDistance(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.RayleighDistance)
      error('You cannot set RayleighDistance property'); 
    end 
   
    function obj = set.k(obj,~)
      fprintf('%s%d\n','k is: ',obj.k)
      error('You cannot set k property'); 
    end 
    
    function obj = set.Waist(obj,~)
      fprintf('%s%d\n','waist is: ',obj.Waist)
      error('You cannot set Waist property'); 
    end 
        
    function obj = set.Radius(obj,~)
      fprintf('%s%d\n','radius is: ',obj.Radius)
      error('You cannot set Waist property'); 
    end 
       
    function obj = set.GouyPhase(obj,~)
      fprintf('%s%d\n','phase is: ',obj.PhiPhase)
      error('You cannot set PhiPhase property'); 
    end 
    
    function obj = set.Amplitude(obj,~)
      fprintf('%s%d\n','amplitude is: ',obj.Amplitude)
      error('You cannot set Amplitude property'); 

    end 
    
    
    %% Get dependent properties
    
    function k = get.k(obj)
       k  = (2*pi)/obj.Wavelength;    
    end
    
    function rayleighDistance = get.RayleighDistance(obj) 
      rayleighDistance  = ((obj.InitialWaist).^2).*(obj.k)/2;    
    end
    
    function Waist = get.Waist(obj)
      Waist = GaussianParameters.waistFunction(obj.PropagationDistance,...
                                               obj.InitialWaist,...
                                               obj.RayleighDistance);
    end
    
    function Phase = get.GouyPhase(obj)
      Phase = GaussianParameters.phaseFunction(obj.PropagationDistance,...
                                               obj.RayleighDistance);
    end
    
    function Radius = get.Radius(obj)
      Radius = GaussianParameters.radiusFunction(obj.PropagationDistance,...
                                                 obj.RayleighDistance);
    end
    
    function Amplitude = get.Amplitude(obj)
      Amplitude = 1./obj.Waist;
    end
    
    function DivergenceAngle = get.DivergenceAngle(obj)
      DivergenceAngle = atan((obj.InitialWaist)./(obj.RayleighDistance));
    end

    
    
  end
  

end