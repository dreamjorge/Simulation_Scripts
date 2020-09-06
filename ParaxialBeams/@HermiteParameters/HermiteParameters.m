classdef HermiteParameters < GaussianParameters & handle & matlab.mixin.Copyable
  
  properties
    n
    m
  end
  
  properties (Dependent)
    HermiteWaistX
    HermiteWaistY
    HermiteWaist
    PhiPhase
  end
  
  methods(Static)
    waistHermite             = getWaist(zDistance,...
                                        InitialWaist,...
                                        RayleighDistance,...
                                        nu,...
                                        mu);
                                      
    waistHermiteOneDirection = getWaistOneDirection(zDistance,...
                                                    InitialWaist,...
                                                    RayleighDistance,...
                                                    l);
  end
      
  methods
    
    function PhiPhase = get.PhiPhase(obj) 
      PhiPhase = (abs(obj.n)+2*obj.m).*obj.GouyPhase;
    end

    function HermiteWaistX = get.HermiteWaistX(obj)
      HermiteWaistX = HermiteParameters.getWaistOneDirection(obj.zCoordinate,obj.InitialWaist,obj.RayleighDistance,obj.n);
    end
   
    function HermiteWaistY = get.HermiteWaistY(obj)
      HermiteWaistY = HermiteParameters.getWaistOneDirection(obj.zCoordinate,obj.InitialWaist,obj.RayleighDistance,obj.m);
    end
        
    function HermiteWaist = get.HermiteWaist(obj)
      HermiteWaist = HermiteParameters.getWaist(obj.zCoordinate,obj.InitialWaist,obj.RayleighDistance,obj.n,obj.m);
    end
    
    %% Constructor
    function Parameters = HermiteParameters(zCoordinate,InitialWaist,Wavelength,nu,mu)
     
     Parameters@GaussianParameters(zCoordinate,InitialWaist,Wavelength);
     
     if nargin == 5 
       Parameters.n             = nu;
       Parameters.m             = mu;
     else
        error('You need introduce PropagationDistance, InitialWaist and Wavelength, l, p inputs')
     end
     
    end
    
  end
end