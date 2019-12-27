classdef HankelLaguerre... < LaguerreBeam & XLaguerreBeam
  
  properties   
    HankelType
    OpticalField
  end
  
%   properties (Dependent)
%     LaguerreWaist
%   end
% %   
  methods(Static)
    [ray] = getLaguerreSlopes(ray,x,y,z,...
                              dx,dy,dz,...
                              xi,yi,zi,...
                              InitialWaist,Wavelength,nu,mu,nh)
  end
%   
  methods
    
%     function LaguerreWaist = get.LaguerreWaist(obj)
%       LaguerreWaist = obj.l+obj.p;
%     end
    
    function Hankel = HankelLaguerre(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu,nh)
    
      
      Hankel.HankelType = nh;
      
      LB  = LaguerreBeam(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu);
      XLB = XLaguerreBeam(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu);
      
      if nh == 1
        Hankel.OpticalField = LB.OpticalField + 1i*XLB.OpticalField;
      elseif nh == 2
        Hankel.OpticalField = LB.OpticalField - 1i*XLB.OpticalField;
      end
      
    end
  end
  
  
end