classdef eLaguerreBeam <  matlab.mixin.Copyable & handle & LaguerreParameters & GaussianBeam
  %% Laguerre Gaussian Beam is a scalar optical field with its parameters defined in
  %properties.
  
  % Example:
  % LB = LaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters);
  % where PropagationDistance,rCoordinate,rCoordinate can be scalar, vector or matrix.  

  properties (Dependent)
    %% Properties dependient of Laguerre Gauss Beam
    LaguerreAmplitude     % Laguerre Gauss factor for Amplitude
    OpticalFieldLaguerre  % Optical Field of Laguerre Gauss
  end
  
  properties
    %% Properties independent of Laguerre Gauss Beam
    thetaCoordinate       % Azimutal coordinate in Cylindrical Coordinates
    Normalization         % Normalization Factor 
  end
  
  methods
    
    function LaguerreAmplitude = get.LaguerreAmplitude(obj)
      %% Factor of Amplitude for Laguerre Gauss Beam
      LaguerreAmplitude = 1;...(1./obj.Waist).*((sqrt(2)*(obj.RadialCoordinate))./obj.Waist).^abs(obj.p);%obj.l);
    end

    function Normalization = get.Normalization(obj)
     %% Factor of Normalizartion for Laguerre Gauss Beam
     
      q = obj.zCoordinate-1i*1i*obj.RayleighDistance;
      Normalization = (-1i*obj.RayleighDistance./q).^(abs(obj.l)+obj.p+1);
    end

    function Laguerre = eLaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters)
      %% Constructor of Laguerre Gauss Beam
      
      % Copying LaguerreParameters to LaguerreBeam Object
      Laguerre@LaguerreParameters(LaguerreParameters.zCoordinate...
                                 ,LaguerreParameters.InitialWaist...
                                 ,LaguerreParameters.Wavelength...
                                 ,LaguerreParameters.l...
                                 ,LaguerreParameters.p);
      
      %Copying Gaussian Beam to LaguerreBeam Object
      Laguerre@GaussianBeam(rCoordinate,LaguerreParameters); 
      % Copying input to properties 
      Laguerre.thetaCoordinate  = thetaCoordinate;

    end
    
    function opticalField =get.OpticalFieldLaguerre(obj)
      %% Calculating Optical Field of Laguerre
      PhiPhase     = obj.PhiPhase;
      l            = obj.l;
      p            = obj.p;
      theta        = obj.thetaCoordinate;
      r            = obj.rCoordinate;
%       xArgument    = (2*(r.^2))./(waist.^2);
      
      q            = obj.zCoordinate+1i*obj.RayleighDistance;
      k            = obj.k;
      alpha        = 1i*(k)./(2*q); 
      
      opticalField = obj.Normalization.*...
                     ...Laguerre.LaguerreAmplitude.*... 
                     exp( 1i*PhiPhase).*...
                     exp(-1i*abs(p)*(theta)).*...
                     LaguerreParameters.getAssociatedLaguerrePolynomial(l,abs(p),alpha.*r.^2).*...
                     obj.OpticalField;
    end
    
    function [] =set.OpticalFieldLaguerre(~,~)
      %% Block OpticalFieldLaguerre as input in Lagurre Gauss Beam
       error('You cannot set OpticalField property'); 
    end
  end
  
end