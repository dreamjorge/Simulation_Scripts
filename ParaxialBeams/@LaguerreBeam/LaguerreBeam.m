classdef LaguerreBeam <  GaussianBeam
  
  properties
    l
    p
  end
  
  properties (Dependent)
    LaguerreWaist
    PhiPhase
    LaguerreAmplitude
  end
  
  properties (Hidden)
    RadialCoordinate
    ThetaCoordinate
    Normalization
  end
  
  methods(Static)
    Lg     = AssociatedLaguerrePolynomial(nu,mu,x);
    waistL = waistLaguerre(zDistance,InitialWaist,RayleighDistance,nu,mu);
  end
  
  methods
    
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (abs(obj.l)+2*obj.p+1).*obj.GouyPhase;
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = LaguerreBeam.waistLaguerre(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      LaguerreAmplitude = (sqrt(2)*(obj.RadialCoordinate)./obj.Waist).^(obj.l);
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function Laguerre = LaguerreBeam(x,y,PropagationDistance, InitialWaist,Wavelength,p,l)
      
      Laguerre@GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength); 
      
      Laguerre.l = l;
      Laguerre.p = p;
      
      [Laguerre.ThetaCoordinate,Laguerre.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      Laguerre.OpticalField = Laguerre.Normalization.*...
                              Laguerre.LaguerreAmplitude.*... 
                              exp(-1i*Laguerre.PhiPhase).*exp(1i*abs(l)*Laguerre.ThetaCoordinate).*...
                              LaguerreBeam.AssociatedLaguerrePolynomial(p,abs(l),(2*Laguerre.RadialCoordinate.^2)./Laguerre.Waist.^2).*...
                              Laguerre.OpticalField;
    end
  end
  
end