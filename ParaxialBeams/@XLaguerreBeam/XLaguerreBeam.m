classdef XLaguerreBeam < GaussianBeam & LaguerreParameters

  properties (Dependent)
    XLaguerreAmplitude
  end
  
  properties (Hidden)
%     RadialCoordinate
%     ThetaCoordinate
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

    
    function XLaguerre = XLaguerreBeam(rCoordinate,thetaCoordinate,zCoordinate, InitialWaist,Wavelength,nu,mu)
      
      
      XLaguerre@GaussianBeam(rCoordinate,zCoordinate,InitialWaist,Wavelength); 
      XLaguerre@LaguerreParameters(zCoordinate,InitialWaist,Wavelength,nu,mu);

%       [XLaguerre.ThetaCoordinate,XLaguerre.RadialCoordinate] = cart2pol(x,y);
%       
      %% Optical Field
      XLaguerre.OpticalField =... XLaguerre.Normalization.*...
                              ...XLaguerre.XLaguerreAmplitude.*... 
                              exp(1i*XLaguerre.PhiPhase).*exp(-1i*abs(mu)*thetaCoordinate).*...
                              XLaguerreBeam.XAssociatedLaguerrePolynomial(nu,abs(mu),(2*rCoordinate.^2)./XLaguerre.Waist.^2).*...
                              XLaguerre.OpticalField;
    end
  end
  

  
end