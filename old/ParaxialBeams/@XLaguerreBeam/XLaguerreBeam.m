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
  
 
  methods


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
      waist        = obj.Waist;
      theta        = obj.thetaCoordinate;
      r            = obj.rCoordinate;
      xArgument    = (2*(r.^2))./(waist.^2);
       
      opticalField =... XLaguerre.Normalization.*...
                   ...XLaguerre.XLaguerreAmplitude.*... 
                    exp(1i*PhiPhase).*exp(-1i*abs(p)*theta).*...
                    LaguerreParameters.getXAssociatedLaguerrePolynomial(l,abs(p),xArgument).*...
                    obj.OpticalField;
     end
  end
  

  
end