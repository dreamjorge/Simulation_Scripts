classdef eHermiteBeam < matlab.mixin.Copyable & handle & HermiteParameters & GaussianBeam
  
  properties
    OpticalFieldHermite
    x
    y
  end
  
  properties (Dependent)
    HermiteAmplitude
  end
  
  properties (Hidden)
    Normalization
  end
  
  methods 

    function HermiteAmplitude = get.HermiteAmplitude(obj)
      HermiteAmplitude = 1;
    end
    
    function Normalization        = get.Normalization(obj)
      Normalization =  1;
    end

    function OpticalFieldHermite  = get.OpticalFieldHermite(obj)
    %% Obatining Optical Field of Hermite
      q = ( 1./obj.Radius-1i./(obj.Waist).^2);
%        q = obj.zCoordinate-1i*obj.RayleighDistance;   
      [Hn,~] = HermiteBeam.getHermiteSolutions(obj.n,sqrt(1i*q).*obj.x);

      [Hm,~] = HermiteBeam.getHermiteSolutions(obj.m,sqrt(1i*q).*obj.y);

      OpticalFieldHermite = obj.Normalization.*...
                            obj.HermiteAmplitude.*...
                            exp(1i*obj.PhiPhase).*...
                            Hn.*Hm.*...
                            obj.OpticalField;
    end
      
    function Hermite          = eHermiteBeam(x,y,hermiteParameters)
    %% Constructor of Hermite
      
      % Copying HermtiteParameters to Hermite Beam Object
      Hermite@HermiteParameters( hermiteParameters.zCoordinate...
                               , hermiteParameters.InitialWaist...
                               , hermiteParameters.Wavelength...
                               , hermiteParameters.n...
                               , hermiteParameters.m);
      
      % Gaussian Beam reqs radial coordinate
      [~,r]=cart2pol(x,y);
      
      % Generate Gaussian Beam for Hermite Beam
      Hermite@GaussianBeam(r,hermiteParameters);
      
      % Copying coordinates to object for generate Optical Field of Hermite
      Hermite.x = x;
      Hermite.y = y;
    
    end
  end
  
end