classdef XLaguerreBeam < GaussianBeam & LaguerreParameters

  properties (Dependent)
    XLaguerreAmplitude
  end
  
  properties (Hidden)
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

    
     function XLaguerre = XLaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters)
      
      XLaguerre@LaguerreParameters(LaguerreParameters.zCoordinate...
                                  ,LaguerreParameters.InitialWaist...
                                  ,LaguerreParameters.Wavelength...
                                  ,LaguerreParameters.l...
                                  ,LaguerreParameters.p);
      
      XLaguerre@GaussianBeam(rCoordinate,LaguerreParameters); 

      PhiPhase =  LaguerreParameters.PhiPhase;
      l        =  LaguerreParameters.p;
      p        =  LaguerreParameters.l;
      waist    =  LaguerreParameters.Waist;

%       [XLaguerre.ThetaCoordinate,XLaguerre.RadialCoordinate] = cart2pol(x,y);
%       
      %% Optical Field
      XLaguerre.OpticalField =... XLaguerre.Normalization.*...
                              ...XLaguerre.XLaguerreAmplitude.*... 
                              exp(1i*PhiPhase).*exp(-1i*abs(p)*thetaCoordinate).*...
                              XLaguerreBeam.XAssociatedLaguerrePolynomial(l,abs(p),(2*rCoordinate.^2)./waist.^2).*...
                              XLaguerre.OpticalField;
    end
  end
  

  
end