classdef GaussianParameters <  handle & matlab.mixin.Copyable
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

% Example: Parameters = GaussianParameters(zCoordinate,InitialWaist,Wavelength)

  properties
    %% Independent properties 
    zCoordinate      % z-coordinate (or distance of propagation)
    InitialWaist     % initial wasit of Gaussian Beam
    Wavelength       % wavelength of Gaussian Beam
  end
  
  properties(Access = private)
     %% Don't set Optical Field
    OpticalField
  end
  
  
  properties (Dependent)
    %% Properties dependent of Gaussian Beam
    RayleighDistance % RayleighDistance.
    k                % k-number.
    Waist            % Waist of Gaussian Beam at zCoordinate.
    Radius           % Radius of Gaussian Beam at zCoordinate.
    GouyPhase        % GouyPhace of Gaussian Beam at zCoordinate.
    Amplitude        % Amplitude of Gaussian Beam at zCoordinate.
    DivergenceAngle  % Angle of Divergence of Gaussian Beam
  end

  
  methods(Static)
    %% functions ofr obtain parameters of Gaussian Beam     
    Waist  = waistFunction(zCoordinate,...
                           InitialWaist,...
                           RayleighDistance);

    Phase  = phaseFunction(zCoordinate,...
                           RayleighDistance);

    Radius = radiusFunction(zCoordinate,...
                            RayleighDistance);

  end
  methods 
    
    function output = copyObject(obj, output)
      %% copying object
       C = metaclass(obj);
       P = C.Properties;
       for kk = 1:length(P)
         if ~P{kk}.Dependent
           output.(P{kk}.Name) = obj.(P{kk}.Name);
         end

       end
     end


    function Parameters = GaussianParameters(zCoordinate,InitialWaist,Wavelength)
     %% Constructor of Gaussian Beam   
      if nargin == 3 
        Parameters.zCoordinate         = zCoordinate;
        Parameters.InitialWaist        = InitialWaist;
        Parameters.Wavelength          = Wavelength;
      else
         error('You need introduce zCoordinate, InitialWaist and Wavelength inputs')
      end
    end
    %% Error for dependent properties
    
    function [] = set.RayleighDistance(obj,~)
      fprintf('%s%d\n','RayleighDistance is: ',obj.RayleighDistance)
      error('You cannot set RayleighDistance property'); 
    end 
   
    function [] = set.k(obj,~)
      fprintf('%s%d\n','k is: ',obj.k)
      error('You cannot set k property'); 
    end 
    
    function [] = set.Waist(obj,~)
      fprintf('%s%d\n','waist is: ',obj.Waist)
      error('You cannot set Waist property'); 
    end 
        
    function [] = set.Radius(obj,~)
      fprintf('%s%d\n','radius is: ',obj.Radius)
      error('You cannot set Waist property'); 
    end 
       
    function [] = set.GouyPhase(obj,~)
      fprintf('%s%d\n','phase is: ',obj.PhiPhase)
      error('You cannot set PhiPhase property'); 
    end 
    
    function [] = set.Amplitude(obj,~)
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