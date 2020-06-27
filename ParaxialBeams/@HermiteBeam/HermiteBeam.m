classdef HermiteBeam < matlab.mixin.Copyable & handle & HermiteParameters & GaussianBeam
  
  properties
    OpticalFieldHermite
    x
    y
  end
  
  properties (Dependent)
%     HermiteWaist
%     HermiteWaistX
%     HermiteWaistY
%     PhiPhase
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
%     function PhiPhase         = get.PhiPhase(obj)
%       PhiPhase = (obj.n+obj.m).*obj.GouyPhase;
%     end

 %   function HermiteWaist     = get.HermiteWaist(obj)
  %    HermiteWaist = HermiteBeam.waistHermite(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.n,obj.m);
  %  end
    
    function HermiteAmplitude = get.HermiteAmplitude(obj)
      HermiteAmplitude = 1;
    end
    
    function Normalization        = get.Normalization(obj)
      Normalization =  1;
    end

    function OpticalFieldHermite  = get.OpticalFieldHermite(obj)

      [Hn,~] = HermiteBeam.hermiteSolutions(obj.n,(sqrt(2)./obj.Waist).*obj.x);
      
      [Hm,~] = HermiteBeam.hermiteSolutions(obj.m,(sqrt(2)./obj.Waist).*obj.y);

      OpticalFieldHermite = obj.Normalization.*...
                            obj.HermiteAmplitude.*...
                            exp(1i*obj.PhiPhase).*...
                            Hn.*Hm.*...
                            obj.OpticalField;
    end
      
    function Hermite          = HermiteBeam(x,y,hermiteParameters)
    %% Constructor of Hermite
      
      % Copying HermtiteParameters to Hermite Beam Object
      Hermite@HermiteParameters( hermiteParameters.zCoordinate...
                               , hermiteParameters.InitialWaist...
                               , hermiteParameters.Wavelength...
                               , hermiteParameters.n...
                               , hermiteParameters.m);
      
      
      [~,r]=cart2pol(x,y);
      
      Hermite@GaussianBeam(r,hermiteParameters);
      
      Hermite.x = x;
      Hermite.y = y;
    
    end
  end
  
end