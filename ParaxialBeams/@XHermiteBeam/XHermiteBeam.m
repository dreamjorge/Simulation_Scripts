classdef XHermiteBeam <  GaussianBeam
  
  properties
    n
    m
  end
  
  properties (Dependent)
    XHermiteWaist
    PhiPhase
    XHermiteAmplitude
  end
  
  properties (Hidden)
    Normalization
  end
  
  methods(Static) 
    [HG,NHG] = hermiteSolutions(nu,x);
  end
  
  methods 
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (obj.n+obj.m).*obj.GouyPhase;
    end

    function XHermiteWaist = get.XHermiteWaist(obj)
      XHermiteWaist = XHermiteBeam.waistXHermite(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end
    
    function XHermiteAmplitude = get.XHermiteAmplitude(obj)
      XHermiteAmplitude = 1;
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  1;...sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function XHermite = XHermiteBeam(x,y,PropagationDistance, InitialWaist,Wavelength,nu,mu)
      
      XHermite@GaussianBeam(x,y,PropagationDistance,InitialWaist,Wavelength); 
      
      XHermite.n = nu;
      XHermite.m = mu;
      [Hn,~]    = XHermite.hermiteSolutions(nu,(sqrt(2)/XHermite.Waist)*x);
      [Hm,~]    = XHermite.hermiteSolutions(mu,(sqrt(2)/XHermite.Waist)*y);

      
      %% Optical Field
      XHermite.OpticalField = XHermite.Normalization.*...
                             XHermite.XHermiteAmplitude.*... 
                             exp(-1i*XHermite.PhiPhase).*...
                             Hn.*Hm.*...
                             XHermite.OpticalField;
    end
  end
  
end