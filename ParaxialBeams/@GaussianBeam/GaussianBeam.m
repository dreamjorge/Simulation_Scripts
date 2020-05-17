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
    
     function cpObj = copyElement(obj)
         % Copy super_prop
         cpObj = copyElement@GaussianParameters(obj);
         % Copy sub_prop1 in subclass

      end
    
    function beam = GaussianBeam(r,gaussianParameters)

      %copying gaussian parameters to beam class
      beam@GaussianParameters(gaussianParameters.zCoordinate...
                             ,gaussianParameters.InitialWaist...
                             ,gaussianParameters.Wavelength);
      ...
      if nargin == 2     

      else 
         error('You need introduce r (value,vector or matrix), and gaussianParameters obect as inputs')
      end 
       
       % Gaussian beam
       beam.OpticalField = beam.Amplitude.*exp(-(r.^2)./(beam.Waist.^2))...
                         .*exp(-1i*beam.k*(r.^2)./(2*beam.Radius))...
                         .*exp(1i*beam.GouyPhase)...
                         .*exp(1i*beam.k*beam.zCoordinate);
       
    end
    

  end
end