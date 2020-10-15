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

    function XHermiteAmplitude = get.XHermiteAmplitude(obj)
      XHermiteAmplitude = 1./(obj.Waist);
    end
    
    function Normalization = get.Normalization(obj)
      Normalization =  1;
    end

    function OpticalFieldHermite  = get.OpticalFieldXHermite(obj)
    %% Obatining Optical Field of Hermite  
      [~,NHn] = HermiteParameters.getHermiteSolutions(obj.n,(sqrt(2)./obj.Waist).*obj.x);

      [~,NHm] = HermiteParameters.getHermiteSolutions(obj.m,(sqrt(2)./obj.Waist).*obj.y); 

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