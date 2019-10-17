classdef LaguerreBeam <  GaussianBeam & LaguerreParameters
  

  properties (Dependent)
    LaguerreAmplitude
  end
  
  properties (Hidden)
    RadialCoordinate
    ThetaCoordinate
    Normalization
  end
  
  methods(Static)
    Lg     = AssociatedLaguerrePolynomial(nu,mu,x);
  end
  
  methods
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      LaguerreAmplitude = (1./obj.Waist).*((sqrt(2)*(obj.RadialCoordinate))./obj.Waist).^abs(obj.l);%obj.l);
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function Laguerre = LaguerreBeam(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu)
      
      Laguerre@GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength); 
      Laguerre@LaguerreParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu)

   
      [Laguerre.ThetaCoordinate,Laguerre.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      Laguerre.OpticalField = Laguerre.Normalization.*...
                              Laguerre.LaguerreAmplitude.*... 
                              exp(1i*Laguerre.PhiPhase).*exp(-1i*abs(mu)*Laguerre.ThetaCoordinate).*...
                              LaguerreBeam.AssociatedLaguerrePolynomial(nu,abs(mu),(2*Laguerre.RadialCoordinate.^2)./Laguerre.Waist.^2).*...
                              Laguerre.OpticalField;
    end
  end
  
end