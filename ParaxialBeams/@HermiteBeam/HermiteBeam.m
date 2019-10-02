classdef HermiteBeam <  GaussianBeam
  
  properties
    n
    m
  end
  
  properties (Dependent)
    HermiteWaist
    PhiPhase
    HermiteAmplitude
  end
  
  properties (Hidden)
    Normalization
  end
  
  methods 
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (abs(obj.l)+2*obj.p+1).*obj.GouyPhase;
    end

    function HermiteWaist = get.HermiteWaist(obj)
      HermiteWaist = HermiteBeam.waistHermite(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end
    
    function HermiteAmplitude = get.HermiteAmplitude(obj)
      HermiteAmplitude = (sqrt(2)*(obj.RadialCoordinate)./obj.Waist).^(obj.l);
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function Hermite = HermiteBeam(x,y,PropagationDistance, InitialWaist,Wavelength,p,l)
      
      Hermite@GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength); 
      
      Hermite.l = l;
      Hermite.p = p;
      
      [Hermite.ThetaCoordinate,Hermite.RadialCoordinate] = cart2pol(x,y);
      
      %% Optical Field
      Hermite.OpticalField = Hermite.Normalization.*...
                              Hermite.HermiteAmplitude.*... 
                              exp(-1i*Hermite.PhiPhase).*exp(1i*abs(l)*Hermite.ThetaCoordinate).*...
                              HermiteBeam.AssociatedHermitePolynomial(p,abs(l),(2*Hermite.RadialCoordinate.^2)./Hermite.Waist.^2).*...
                              Hermite.OpticalField;
    end
  end
  
end