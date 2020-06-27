classdef XHermiteBeam < matlab.mixin.Copyable & handle & HermiteParameters & GaussianBeam
  
  properties
    OpticalFieldXHermite
    x
    y
  end
  
  properties (Dependent)
    XHermiteAmplitude
  end
  
  properties (Hidden)
    Normalization
  end
  
  methods(Static) 
    [HG,NHG] = hermiteSolutions(nu,x);
  end
  
  methods 
%     function PhiPhase = get.PhiPhase(obj)
%       PhiPhase = (obj.n+obj.m).*obj.GouyPhase;
%     end
% 
%     function XHermiteWaist = get.XHermiteWaist(obj)
%       XHermiteWaist = XHermiteBeam.waistXHermite(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
%     end
%     
    function XHermiteAmplitude = get.XHermiteAmplitude(obj)
      XHermiteAmplitude = 1;
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  1;...sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function OpticalFieldHermite  = get.OpticalFieldXHermite(obj)
      
      [~,NHn] = HermiteBeam.hermiteSolutions(obj.n,(sqrt(2)./obj.Waist).*obj.x);

      [~,NHm] = HermiteBeam.hermiteSolutions(obj.m,(sqrt(2)./obj.Waist).*obj.y); 

      OpticalFieldHermite = obj.Normalization.*...
                            obj.XHermiteAmplitude.*... 
                            exp(1i*obj.PhiPhase).*...
                            NHn.*NHm.*...
                            obj.OpticalField;
    end  
  
  
    function XHermite = XHermiteBeam(x,y,hermiteParameters)
      
      % Copying HermtiteParameters to Hermite Beam Object
      XHermite@HermiteParameters( hermiteParameters.zCoordinate...
                                , hermiteParameters.InitialWaist...
                                , hermiteParameters.Wavelength...
                                , hermiteParameters.n...
                                , hermiteParameters.m);
      
      
      [~,r]=cart2pol(x,y);
      
      XHermite@GaussianBeam(r,hermiteParameters);
      
      XHermite.x = x;
      XHermite.y = y;
    
    end
    
  end
  
end