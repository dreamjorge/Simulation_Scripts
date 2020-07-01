classdef HermiteParameters < GaussianParameters
  
  properties
    n
    m
  end
  
  properties (Dependent)
    HermiteWaist
    PhiPhase
  end
  
  methods(Static)
    waistL = Waist(zDistance,InitialWaist,RayleighDistance,nu,mu);
  end
    
   methods
    
    function PhiPhase = get.PhiPhase(obj) 
      PhiPhase = (abs(obj.n)+2*obj.m).*obj.GouyPhase;
    end

    function HermiteWaist = get.HermiteWaist(obj)
      HermiteWaist = LaguerreParameters.Waist(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.n,obj.m);
    end
   
    %% Constructor
    function Parameters = HermiteParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu)
     
     Parameters@GaussianParameters(PropagationDistance,InitialWaist,Wavelength);
     
     if nargin == 5 
       Parameters.n             = nu;
       Parameters.m             = mu;
     else
        error('You need introduce PropagationDistance, InitialWaist and Wavelength, l, p inputs')
     end

    end
  end
end