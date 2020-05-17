classdef GaussianParameters <  handle
% This class gives Gaussian parameters of Gaussian Beam
% recives zCoordinate,InitialWaist,Wavelength as input and
% gives a object with next properties:
%
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
% Example: Parameters = GaussianParameters(zCoordinate,InitialWaist,Wavelength)

  properties
    zCoordinate    % z-coordinate (or distance of propagation)
    InitialWaist   % initial wasit of Gaussian Beam
    Wavelength     % wavelength of Gaussian Beam
  end
  
  properties (Dependent)
    RayleighDistance % RayleighDistance.
    k                % k-number.
    Waist            % Waist of Gaussian Beam at zCoordinate.
    Radius           % Radius of Gaussian Beam at zCoordinate.
    GouyPhase        % GouyPhace of Gaussian Beam at zCoordinate.
    Amplitude        % Amplitude of Gaussian Beam at zCoordinate.
    DivergenceAngle  % Angle of Divergence of Gaussian Beam
  end

  
  methods(Static)
%     
        Waist  = waistFunction(zCoordinate,...
                               InitialWaist,...
                               RayleighDistance);
                                    
        Phase  = phaseFunction(zCoordinate,...
                               RayleighDistance);
                                         
        Radius = radiusFunction(zCoordinate,...
                                RayleighDistance);

  end
  
  
  methods
    
  %% Constructor
     function Parameters = GaussianParameters(zCoordinate,InitialWaist,Wavelength)
      if nargin == 3 
        Parameters.zCoordinate         = zCoordinate;
        Parameters.InitialWaist        = InitialWaist;
        Parameters.Wavelength          = Wavelength;
      else
         error('You need introduce zCoordinate, InitialWaist and Wavelength inputs')
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
      Waist = GaussianParameters.waistFunction(obj.zCoordinate,...
                                               obj.InitialWaist,...
                                               obj.RayleighDistance);
    end
    
    function Phase = get.GouyPhase(obj)
      Phase = GaussianParameters.phaseFunction(obj.zCoordinate,...
                                               obj.RayleighDistance);
    end
    
    function Radius = get.Radius(obj)
      Radius = GaussianParameters.radiusFunction(obj.zCoordinate,...
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