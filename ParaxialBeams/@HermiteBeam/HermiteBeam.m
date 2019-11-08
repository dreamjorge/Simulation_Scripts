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
  
  methods(Static) 
    [HG,NHG] = hermiteSolutions(nu,x);
    wh       = waistHermite(zDistance,initialWaist,rayleighDistance,n,m);
  end
  
  methods 
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (obj.n+obj.m).*obj.GouyPhase;
    end

    function HermiteWaist = get.HermiteWaist(obj)
      HermiteWaist = HermiteBeam.waistHermite(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.n,obj.m);
    end
    
    function HermiteAmplitude = get.HermiteAmplitude(obj)
      HermiteAmplitude = 1;
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  1;
    end

    function Hermite = HermiteBeam(x,y,PropagationDistance, InitialWaist,Wavelength,n,m)
      
      Hermite@GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength); 
      
      Hermite.n = n;
      Hermite.m = m;
      [Hn,~]    = Hermite.hermiteSolutions(n,(sqrt(2)./Hermite.Waist).*x);
      [Hm,~]    = Hermite.hermiteSolutions(m,(sqrt(2)./Hermite.Waist).*y);

      
      %% Optical Field
      Hermite.OpticalField = Hermite.Normalization.*...
                             Hermite.HermiteAmplitude.*... 
                             exp(1i*Hermite.PhiPhase).*...
                             Hn.*Hm.*...
                             Hermite.OpticalField;
    end
  end
  
end