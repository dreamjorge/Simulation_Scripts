classdef HankelHermite... < LaguerreBeam & XLaguerreBeam
  
  properties   
    HankelType
    OpticalField
  end
  
  methods(Static)
    [ray] = getHermiteSlopes(ray,x,y,z,...
                             dx,dy,dz,...
                             xi,yi,zi,...
                             InitialWaist,Wavelength,nu,mu,nh)
  end
  
  methods
    
    function Hankel = HankelHermite(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu,xnh,ynh)

      Hankel.HankelType = [xnh,ynh];

      HB  =  HermiteBeam(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu);
      XHB = XHermiteBeam(x,y,PropagationDistance,InitialWaist,Wavelength,nu,mu);

      if     ((xnh == 1) && (ynh == 1))
        Hankel.OpticalField = (HB.OpticalField + 1i*XHB.OpticalField)...
                            .*(HB.OpticalField + 1i*XHB.OpticalField);

      elseif ((xnh == 1) && (ynh == 2))
        Hankel.OpticalField = (HB.OpticalField + 1i*XHB.OpticalField)...
                            .*(HB.OpticalField - 1i*XHB.OpticalField);

      elseif ((xnh == 2) && (ynh == 1))
        Hankel.OpticalField = (HB.OpticalField - 1i*XHB.OpticalField)...
                            .*(HB.OpticalField + 1i*XHB.OpticalField);

      elseif ((xnh == 2) && (ynh ==2))
        Hankel.OpticalField = (HB.OpticalField - 1i*XHB.OpticalField)...
                            .*(HB.OpticalField - 1i*XHB.OpticalField);

      end

    end
    
  end  
  
  
end