classdef HankelHermite
  
  properties   
    HankelType
    OpticalField
  end
  
  methods(Static)
    Rays = getPropagateCartesianRays(Rays,...
                                     x,y,...
                                     difr,...
                                     HParametersZi,...
                                     HParametersZ,...
                                     HankelType) 
  end
  
  methods
    
    function Hankel = HankelHermite(x,y,hermiteParameters,hankelType)

      Hankel.HankelType = hankelType;
      WaistGauss        = hermiteParameters.Waist;

      [Hx,NHx] = ...
      HermiteBeam.hermiteSolutions(hermiteParameters.n...
                                  ,(sqrt(2)./WaistGauss).*x);
      
      [Hy,NHy] = ...
      HermiteBeam.hermiteSolutions(hermiteParameters.m...
                                  ,(sqrt(2)./WaistGauss).*y);
                            
      
      GaussB   = GaussianBeam(sqrt(x.^2+y.^2),hermiteParameters).OpticalField;


      if     (hankelType == 11)
        Hankel.OpticalField = (Hx+1i*NHx).*(Hy+1i*NHy).*GaussB;
       
      elseif (hankelType == 12)
        Hankel.OpticalField = (Hx+1i*NHx).*(Hy-1i*NHy).*GaussB;

      elseif (hankelType == 21)
        Hankel.OpticalField = (Hx-1i*NHx).*(Hy+1i*NHy).*GaussB;

      elseif (hankelType == 22)
        Hankel.OpticalField = (Hx-1i*NHx).*(Hy-1i*NHy).*GaussB;

      end

    end
    
  end  
  
  
end