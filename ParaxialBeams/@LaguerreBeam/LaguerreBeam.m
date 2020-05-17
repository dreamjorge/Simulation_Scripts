classdef LaguerreBeam < GaussianBeam & LaguerreParameters
  % Laguerre Gaussian Beam is a scalar optical field with its parameters defined in
  % properties.
  % Example:
  % LB = LaguerreBeam(X,Y,PropagationDistance,RayleighDistance,InitialWaist,l,p);
  % where PropagationDistance,X,Y can be scalar, vector or matrix.  

  properties (Dependent)
    LaguerreAmplitude
  end
  
  properties (Hidden)
%     RadialCoordinate
%     ThetaCoordinate
    Normalization
  end
  
  methods(Static)
    Lg     = AssociatedLaguerrePolynomial(nu,mu,x);
  end
  
  methods
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      LaguerreAmplitude = 1;...(1./obj.Waist).*((sqrt(2)*(obj.RadialCoordinate))./obj.Waist).^abs(obj.p);%obj.l);
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function Laguerre = LaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters)
      
      Laguerre@LaguerreParameters(LaguerreParameters.zCoordinate...
                                 ,LaguerreParameters.InitialWaist...
                                 ,LaguerreParameters.Wavelength...
                                 ,LaguerreParameters.l...
                                 ,LaguerreParameters.p);
      
      Laguerre@GaussianBeam(rCoordinate,LaguerreParameters); 

      %% Optical Field
      PhiPhase =  LaguerreParameters.PhiPhase;
      l       =  LaguerreParameters.p;
      p       =  LaguerreParameters.l;
      waist    =  LaguerreParameters.Waist;
      
      Laguerre.OpticalField = ...Laguerre.Normalization.*...
                              ...Laguerre.LaguerreAmplitude.*... 
                              exp(1i*PhiPhase).*exp(-1i*abs(p)*thetaCoordinate).*...
                              LaguerreBeam.AssociatedLaguerrePolynomial(l,abs(p),(2*rCoordinate.^2)./(waist.^2)).*...
                              Laguerre.OpticalField;
    end
  end
  
end