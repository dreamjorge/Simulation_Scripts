classdef LaguerreParameters < GaussianParameters
  
  properties
    l
    p
  end
  
  properties (Dependent)
    LaguerreWaist
    PhiPhase
  end
  
  methods(Static)
    waistL = Waist(zDistance,InitialWaist,RayleighDistance,nu,mu);
    GP = GaussianParameters(PropagationDistance,InitialWaist,Wavelength);
  end
    
   methods
    
    function PhiPhase = get.PhiPhase(obj)
      PhiPhase = (abs(obj.l)+2*obj.p).*obj.GouyPhase;
    end

    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = LaguerreParameters.Waist(obj.PropagationDistance,obj.InitialWaist,obj.RayleighDistance,obj.l,obj.p);
    end
   
    %% Constructor
    function Parameters = LaguerreParameters(PropagationDistance,InitialWaist,Wavelength,nu,mu)
     
     Parameters@GaussianParameters(PropagationDistance,InitialWaist,Wavelength);
     
     if nargin == 5 
       Parameters.l             = nu;
       Parameters.p             = mu;
       %Parameters.LaguerreWaist = get.LaguerreWaist(PropagationDistance,InitialWaist,Parameters.RayleighDistance,nu,mu);
       %Parameters.PhiPhase      = get.PhiPhase(l,p,Parameters.GouyPhase);
     else
        error('You need introduce PropagationDistance, InitialWaist and Wavelength, l, p inputs')
     end

    end
  end
end