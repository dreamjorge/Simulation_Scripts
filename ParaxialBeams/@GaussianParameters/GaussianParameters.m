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
    %% functions for obtain parameters of Gaussian Beam
    Waist  = getWaist(zCoordinate,...
                      InitialWaist,...
                      RayleighDistance);

    Phase  = getPhase(zCoordinate,...
                      RayleighDistance);

    Radius = getRadius(zCoordinate,...
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
        % Input validation
        validateattributes(zCoordinate, {'numeric'}, {'scalar'}, 'GaussianParameters', 'zCoordinate');
        validateattributes(InitialWaist, {'numeric'}, {'positive', 'scalar'}, 'GaussianParameters', 'InitialWaist');
        validateattributes(Wavelength, {'numeric'}, {'positive', 'scalar'}, 'GaussianParameters', 'Wavelength');
        
        Parameters.zCoordinate         = zCoordinate;
        Parameters.InitialWaist        = InitialWaist;
        Parameters.Wavelength          = Wavelength;
      else
         error('You need introduce zCoordinate, InitialWaist and Wavelength inputs')
      end
    end
    
    function str = toString(obj)
      %TOSTRING Return string representation of parameters
      %   str = toString() returns formatted string with all parameters
      %
      %   Example:
      %       params = GaussianParameters(0, 100e-6, 632.8e-9);
      %       disp(toString(params));
      
      str = sprintf([...
          'GaussianParameters:\n', ...
          '  zCoordinate:      %g\n', ...
          '  InitialWaist:     %g\n', ...
          '  Wavelength:       %g\n', ...
          '  RayleighDistance: %g\n', ...
          '  k:                %g\n', ...
          '  Waist:            %g\n', ...
          '  Radius:           %g\n', ...
          '  GouyPhase:        %g\n', ...
          '  DivergenceAngle:  %g\n'], ...
          obj.zCoordinate, obj.InitialWaist, obj.Wavelength, ...
          obj.RayleighDistance, obj.k, obj.Waist, ...
          obj.Radius, obj.GouyPhase, obj.DivergenceAngle);
    end
    
    function bool = isEqual(obj, other)
      %ISEQUAL Compare two GaussianParameters objects
      %   bool = isEqual(other) returns true if all properties match
      
      if ~isa(other, 'GaussianParameters')
        bool = false;
        return;
      end
      
      tolerance = 1e-12;
      bool = abs(obj.zCoordinate - other.zCoordinate) < tolerance && ...
             abs(obj.InitialWaist - other.InitialWaist) < tolerance && ...
             abs(obj.Wavelength - other.Wavelength) < tolerance;
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
      Waist = GaussianParameters.getWaist(obj.zCoordinate,...
                                          obj.InitialWaist,...
                                          obj.RayleighDistance);
    end

    function Phase = get.GouyPhase(obj)
      Phase = GaussianParameters.getPhase(obj.zCoordinate,...
                                          obj.RayleighDistance);
    end

    function Radius = get.Radius(obj)
      Radius = GaussianParameters.getRadius(obj.zCoordinate,...
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