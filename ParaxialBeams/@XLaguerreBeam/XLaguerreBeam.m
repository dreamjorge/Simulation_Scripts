classdef XLaguerreBeam < GaussianBeam & LaguerreParameters

  properties (Dependent)
    XLaguerreAmplitude
  end
  
  properties (Hidden)
    RadialCoordinate
    ThetaCoordinate
    Normalization
  end
  
  methods(Static)
    XLg     = XAssociatedLaguerrePolynomial(nu,mu,x);
  end
  
  methods
    
    
    function XLaguerreAmplitude = get.XLaguerreAmplitude(obj)
      XLaguerreAmplitude = 1;...(sqrt(2)*(obj.RadialCoordinate)./obj.Waist).^(obj.l);
    end

    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    
    function XLaguerre = XLaguerreBeam(x,y,PropagationDistance, InitialWaist,Wavelength,nu,mu)
      
      
      XLaguerre@GaussianBeam(x,y,PropagationDistance, InitialWaist,Wavelength); 
      XLaguerre@LaguerreParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu);

      [XLaguerre.ThetaCoordinate,XLaguerre.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      XLaguerre.OpticalField =... XLaguerre.Normalization.*...
                              ...XLaguerre.XLaguerreAmplitude.*... 
                              exp(1i*XLaguerre.PhiPhase).*exp(-1i*abs(mu)*XLaguerre.ThetaCoordinate).*...
                              XLaguerreBeam.XAssociatedLaguerrePolynomial(nu,abs(mu),(2*XLaguerre.RadialCoordinate.^2)./XLaguerre.Waist.^2).*...
                              XLaguerre.OpticalField;
    end
  end
  
  
  
  
  
end