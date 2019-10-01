classdef GaussianBeam
  % Gaussian Beam is a scalar optical field with its parameters defined in
  % properties.
  % Example:
  % GB = GaussianBeam(X,Y,PropagationDistance,RayleighDistance,InitialWaist);
  % where PropagationDistance,X,Y can be scalar, vector or matrix.
  properties
    GuassianParameters = struct('RayleighDistance',{}...
                               ,'InitialWaist',{}...
                               ,'PropagationDistance',{}...
                               ,'Wavelength',{}...
                               ,'k',{}...
                               ,'Waist',{}...
                               ,'Radius',{}...
                               ,'Phase',{}...
                               ,'Amplitude',{}...
                               );
    OpticalField;
  end
  
  methods
    
    function beam = GaussianBeam(x,y,PropagationDistance,RayleighDistance,InitialWaist)

      ...
      if nargin == 5     
      
        %Calculate parameters of beam
        k            = 2*RayleighDistance/InitialWaist;
        WaveLength   = 2*pi/k;
        Waist        = waistPhysicalGaussianBeam(PropagationDistance,InitialWaist,RayleighDistance);
        Radius       = radiusPhysicalGaussianBeam(PropagationDistance,RayleighDistance);
        Phase        = phasePhysicalGaussianBeam(PropagationDistance,RayleighDistance);
        OpticalField =  (InitialWaist./Waist).*exp(-(x.^2+y.^2)./(Waist).^2)...
          .*exp(k*(x.^2+y.^2)/(2*(Radius))).*exp(-1i*k*PropagationDistance).*exp(-Phase.*PropagationDistance);
        
        %Assign Values to properties
        beam.GuassianParameters(1).InitialWaist        = InitialWaist;
        beam.GuassianParameters(1).RayleighDistance    = RayleighDistance;
        beam.GuassianParameters(1).PropagationDistance = PropagationDistance;
        beam.GuassianParameters(1).Wavelength          = WaveLength;
        beam.GuassianParameters(1).k                   = k;
        beam.GuassianParameters(1).Waist               = Waist;
        beam.GuassianParameters(1).Radius              = Radius;
        beam.GuassianParameters(1).Phase               = Phase;
        beam.OpticalField                              = OpticalField;
      end 
    end
    

  end
end