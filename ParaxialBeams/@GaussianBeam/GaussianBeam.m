classdef GaussianBeam < GaussianParameters
  % Gaussian Beam is a scalar optical field with its parameters defined in
  % properties.
  % Example:
  % GB = GaussianBeam(X,Y,PropagationDistance,RayleighDistance,InitialWaist);
  % where PropagationDistance,X,Y can be scalar, vector or matrix.
  properties
    OpticalField;
  end
  
  methods
    
    function beam = GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength)

      ...
      if nargin == 5     
         super_args{1} = PropagationDistance;
         super_args{2} = InitialWaist;
         super_args{3} = Wavelength;
      else
         error('You need introduce x, y, PropagationDistance, InitialWaist and Wavelength inputs')
      end 
       beam@GaussianParameters(super_args{:});
       
       beam.OpticalField = beam.Amplitude.*exp(-(x.^2+y.^2)./(beam.Waist.^2))...
                         .*exp(-1i*beam.k*(x.^2+y.^2)./(2*beam.Radius))...
                         .*exp(1i*beam.GouyPhase)...
                         .*exp(1i*beam.k*PropagationDistance);
       
    end
    

  end
end