classdef HankelLaguerre 
  
  properties   
    OpticalField
    HankelType
    l
    p
  end
  
  properties (Dependent)
    LaguerreWaist
  end
  
  methods(Static)
    [ray] = getLaguerreSlopes(ray,x,y,z,...
                              dx,dy,dz,...
                              xi,yi,zi,...
                              InitialWaist,Wavelength,p,l,nh)
  end
  
  methods
    
    function LaguerreWaist = get.LaguerreWaist(obj)
      LaguerreWaist = obj.l+obj.p;

    end
    
    function Hankel = HankelLaguerre(x,y,PropagationDistance,InitialWaist,Wavelength,l,p,nh)
    
      
      Hankel.l          = l;
      Hankel.p          = p;
      Hankel.HankelType = nh;
      
      LB  = LaguerreBeam(x,y,PropagationDistance,InitialWaist,Wavelength,l,p);
      XLB = XLaguerreBeam(x,y,PropagationDistance,InitialWaist,Wavelength,l,p);
      
      if nh == 1
        Hankel.OpticalField = LB.OpticalField + 1i*XLB.OpticalField;
      elseif nh == 2
        Hankel.OpticalField = LB.OpticalField - 1i*XLB.OpticalField;
      end
      
    end
  end
  
  
end