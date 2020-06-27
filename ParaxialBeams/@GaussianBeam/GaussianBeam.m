classdef GaussianBeam < matlab.mixin.Copyable & GaussianParameters
  % Gaussian Beam is a scalar optical field with its parameters defined in next
  % properties:
  %  
  % -OpticalField; 
  % -rCoordinate;
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
  
  %% Example:
  % GB = GaussianBeam(r,gaussianParameters);
  % where r can be scalar, vector or matrix; and gaussianParameters defined in 
  % @GaussianParameters
  %%
  properties (Dependent)
    OpticalField; % Optical Field of Gaussian Beam
  end
  
  properties
    rCoordinate;  % Radial Coordinates
  end
  
  methods (Access = protected)

  end
  
  methods 

    function output = copyObject(obj, output)
    %% Copy Object Function
       C = metaclass(obj);
       P = C.Properties;
       for k = 1:length(P)
         if ~P{k}.Dependent
           output.(P{k}.Name) = obj.(P{k}.Name);
         end
       end
    end

    function opticalField = get.OpticalField(obj)
    %% Funtion for get Optical Field of Gaussian Beams
      opticalField = obj.Amplitude...
                  .* exp(-   ((obj.rCoordinate).^2)./((obj.Waist).^2))...
                  .* exp(-1i*((obj.k)*((obj.rCoordinate).^2)./(2*(obj.Radius))))...
                  .* exp( 1i*((obj.GouyPhase)))...
                  .* exp( 1i*((obj.k)*(obj.zCoordinate)));
    end   
    
    function beam = GaussianBeam(rCoordinate,gaussianParameters)
    %% Constructor of Gaussian Beam
      if nargin == 2     
%         if  isa(gaussianParameters,'GaussianParameters')
%           % Do nothing, constructor can be done
%         else
%           error('second input is not a GaussianParameters object')
%         end

      else 
        error('You need introduce r (value,vector or matrix), and gaussianParameters obect as inputs')
      end 
          %Copying Properties of Gaussian Parameters to GaussianBeam
          beam@GaussianParameters(gaussianParameters.zCoordinate...
                                 ,gaussianParameters.InitialWaist...
                                 ,gaussianParameters.Wavelength);
          %Copying Radial coordinate of Input
          beam.rCoordinate = rCoordinate;  
          
          if  isa(gaussianParameters,'LaguerreParameters')
%              beam.OpticalFieldLaguerre = beam.OpticalField;
          end
    end
   
    function [] = set.OpticalField(~,~)
      %fprintf('%s%d\n','OpticalField is: ',obj.OpticalField)
      error('You cannot set OpticalField property'); 
    end   
    
  end
  
end