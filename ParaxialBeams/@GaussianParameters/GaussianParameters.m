classdef GaussianParameters <  handle
  
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
  end

  
  methods(Static)
%     
        waist = waistPhysicalGaussianBeam(PropagationDistance,...
                                      InitialWaist,...
                                      RayleighDistance);
                                    
        phase =  phasePhysicalGaussianBeam(PropagationDistance,...
                                           RayleighDistance);
                                         
        radius = radiusPhysicalGaussianBeam(PropagationDistance,...
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
      fprintf('%s%d\n','RayleighDistance is: ',obj.k)
      error('You cannot set k property'); 
    end 
    
    function obj = set.Waist(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.Waist)
      error('You cannot set Waist property'); 
    end 
        
    function obj = set.Radius(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.Radius)
      error('You cannot set Waist property'); 
    end 
       
    function obj = set.GouyPhase(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.PhiPhase)
      error('You cannot set PhiPhase property'); 
    end 
    
    function obj = set.Amplitude(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.Amplitude)
      error('You cannot set Amplitude property'); 

    end 
    
    
    %% Get dependent properties
    
    function k = get.k(obj)
      
       k  = (2*pi)/obj.Wavelength;    

    end
    
    function rayleighDistance = get.RayleighDistance(obj) 
      rayleighDistance  = obj.InitialWaist.*(obj.k)/2;    
    end
    
    function waist = get.Waist(obj)
      waist = GaussianParameters.waistPhysicalGaussianBeam(obj.PropagationDistance,...
                                        obj.InitialWaist,...
                                        obj.RayleighDistance);
    end
    
    function Phase = get.GouyPhase(obj)
      Phase = GaussianParameters.phasePhysicalGaussianBeam(obj.PropagationDistance,...
                                        obj.RayleighDistance);
    end
    
    function Radius = get.Radius(obj)
      Radius = GaussianParameters.radiusPhysicalGaussianBeam(obj.PropagationDistance,...
                                          obj.RayleighDistance);
    end
    
    function Amplitude = get.Amplitude(obj)
      Amplitude = 1./obj.Waist;
    end
    

    
    
  end
  

end