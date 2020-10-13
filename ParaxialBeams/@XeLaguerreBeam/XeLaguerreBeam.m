classdef XLaguerreBeam < GaussianBeam & LaguerreParameters
  
  properties
    thetaCoordinate
  end
  
  properties (Dependent)
    XLaguerreAmplitude
    OpticalFieldLaguerre
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

      XLaguerre.thetaCoordinate = thetaCoordinate;

%       [XLaguerre.ThetaCoordinate,XLaguerre.RadialCoordinate] = cart2pol(x,y);
%       
      %% Optical Field

     end
    
     
    function [opticalField]=get.OpticalFieldLaguerre(obj)
      PhiPhase     = obj.PhiPhase;
      l            = obj.l;
      p            = obj.p;
      theta        = obj.thetaCoordinate;
      r            = obj.rCoordinate;

      q            = obj.zCoordinate+1i*obj.RayleighDistance;
      k            = obj.k;
      alpha        = 1i*(k)./(2*q); 
      
      
      opticalField =... XLaguerre.Normalization.*...
                   ...XLaguerre.XLaguerreAmplitude.*... 
                    exp(1i*PhiPhase).*exp(-1i*abs(p)*theta).*...
                    XLaguerreBeam.XAssociatedLaguerrePolynomial(l,abs(p),alpha.*r.^2).*...
                    obj.OpticalField;
    end
  end
  

  
end