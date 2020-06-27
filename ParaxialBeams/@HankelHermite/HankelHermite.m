classdef HankelHermite
  
  properties   
    HankelType
    OpticalField
  end
  
  methods(Static)
    Rays = getPropagateCartesianRays(Rays,...
                                     TotalRays,...
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
      PhiPhase          = hermiteParameters.PhiPhase;

      [Hx,NHx] = ...
      HermiteBeam.hermiteSolutions(hermiteParameters.n...
                                  ,(sqrt(2)./WaistGauss).*x);
      
      [Hy,NHy] = ...
      HermiteBeam.hermiteSolutions(hermiteParameters.m...
                                  ,(sqrt(2)./WaistGauss).*y);
      GaussX   = GaussianBeam(x,hermiteParameters).OpticalField;
      GaussY   = GaussianBeam(y,hermiteParameters).OpticalField;
      
      Hx       = Hx .*GaussX.*exp(1i*PhiPhase);
      Hy       = Hy .*GaussY.*exp(1i*PhiPhase);
      NHx      = NHx.*GaussX.*exp(1i*PhiPhase);
      NHy      = NHy.*GaussY.*exp(1i*PhiPhase);

      if     (hankelType == 11)
        Hankel.OpticalField = (Hx+1i*NHx).*(Hy+1i*NHy);
       
      elseif (hankelType == 12)
        Hankel.OpticalField = (Hx+1i*NHx).*(Hy-1i*NHy);

      elseif (hankelType == 21)
        Hankel.OpticalField = (Hx-1i*NHx).*(Hy+1i*NHy);

      elseif (hankelType == 22)
        Hankel.OpticalField = (Hx-1i*NHx).*(Hy-1i*NHy);

      end

    end
    
  end  
  
  
end