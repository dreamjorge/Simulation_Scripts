classdef XeHermiteBeam < matlab.mixin.Copyable & handle & HermiteParameters & GaussianBeam
  
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

    function XHermiteAmplitude = get.XHermiteAmplitude(obj)
      XHermiteAmplitude = 1;
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  1;...sqrt(2*factorial(obj.p)/(pi*factorial(obj.p+abs(obj.l))));
    end

    function OpticalFieldHermite  = get.OpticalFieldXHermite(obj)
    %% Obatining Optical Field of Hermite  
      q = ( 1./obj.Radius-1i./(obj.Waist).^2);
    
      [~,NHn] = HermiteParameters.getHermiteSolutions(obj.n,sqrt(1i*q).*obj.x);

      [~,NHm] = HermiteParameters.getHermiteSolutions(obj.m,sqrt(1i*q).*obj.y); 

      OpticalFieldHermite = obj.Normalization.*...
                            obj.XHermiteAmplitude.*... 
                            exp(1i*obj.PhiPhase).*...
                            NHn.*NHm.*...
                            obj.OpticalField;
    end  
  
  
    function XHermite = XeHermiteBeam(x,y,hermiteParameters)
      
      % Copying HermtiteParameters to Hermite Beam Object
      XHermite@HermiteParameters( hermiteParameters.zCoordinate...
                                , hermiteParameters.InitialWaist...
                                , hermiteParameters.Wavelength...
                                , hermiteParameters.n...
                                , hermiteParameters.m);
                              
      % Gaussian Beam reqs radial coordinate
      [~,r]=cart2pol(x,y);
      
      % Generate Gaussian Beam for Hermite Beam
      XHermite@GaussianBeam(r,hermiteParameters);
      
      % Copying coordinates to object for generate Optical Field of Hermite
      XHermite.x = x;
      XHermite.y = y;
    
    end
    
  end
  
end