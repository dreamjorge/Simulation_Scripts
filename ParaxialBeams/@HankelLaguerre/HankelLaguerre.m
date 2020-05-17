classdef HankelLaguerre 
  
  properties   
    HankelType    % hankel type (1,2)
    OpticalField  % Optical field
  end
    
  methods(Static)
    %Function for obtain slopes in a point given
    [ray] = getLaguerreSlopes(ray,...
                              x,y,z,...
                              dr,...
                              InitialWaist,Wavelength,nu,mu,nh)
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
        Hankel.OpticalField = LB.OpticalField + 1i*XLB.OpticalField;
      elseif nh == 2
        Hankel.OpticalField = LB.OpticalField - 1i*XLB.OpticalField;
      end
      
    end
    
  end
  
  
end