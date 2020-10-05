classdef HankelLaguerre 
  
  properties   
    HankelType    % hankel type (1,2)
    OpticalFieldLaguerre
  end
    
  methods(Static)
    %Function for obtain slopes in a point given
    [ray] = getPropagateCylindricalRays(Rays,...
                                        TotalRays,...
                                        r,th,...
                                        difr,...
                                        LParametersZi,...
                                        LParametersZ,...
                                        HankelType) 
  end
   
  methods
    
    % Constructor of Object
    function Hankel = HankelLaguerre(rCoordinate,thetaCoordinate,LaguerreParameters,nh)
     
      %Hankel Type of input
      Hankel.HankelType = nh;
      %First and second solutions of Laguerre
      LB  =  LaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters);
      XLB = XLaguerreBeam(rCoordinate,thetaCoordinate,LaguerreParameters);
       
      %Taking different sign according type of Hankel
      if nh == 1
        Hankel.OpticalFieldLaguerre = LB.OpticalFieldLaguerre + 1i*XLB.OpticalFieldLaguerre;
      elseif nh == 2
        Hankel.OpticalFieldLaguerre = LB.OpticalFieldLaguerre - 1i*XLB.OpticalFieldLaguerre;
      end
      
    end
    
  end
  
  
end