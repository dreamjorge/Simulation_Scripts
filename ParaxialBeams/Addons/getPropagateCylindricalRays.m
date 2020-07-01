function [Rays] = getPropagateCylindricalRays(Rays,...
                                              TotalRays,...
                                              r,th,...
                                              difr,...
                                              LParametersZi,...
                                              LParametersZ,...
                                              HankelType) 
          
%% Obtain New coordinates of ray using a step of ray in z-direction,
%equation or ray r(z) = mrz*dz+r(z-1) using previous value
debugMode = 'OFF';                                                         
  % each component of diferential dr 
  tempdr     = num2cell(difr);
  [dr,dth,dz] = deal(tempdr{:});
 
  for point_index = 1 : TotalRays

    temporalRay = CylindricalRay(); 

    % copy data a temporal ray 
    temporalRay = copyArrayRay2Ray(Rays,temporalRay,point_index);
    % new coordinates of ray
    temporalRay.rCoordinate     = temporalRay.rCoordinate     + (1./temporalRay.zrSlope )*dz;
    temporalRay.thetaCoordinate = temporalRay.thetaCoordinate + (1./temporalRay.zthSlope)*dth;
    temporalRay.zCoordinate     = temporalRay.zCoordinate + dz;
    % obtain new (x,y) coordinates 
    [temporalRay.xCoordinate,temporalRay.yCoordinate] = pol2cart(temporalRay.thetaCoordinate,temporalRay.rCoordinate);

    % taking this new coordinates of point to calculate hankels
    ri  = temporalRay.rCoordinate;
    thi = temporalRay.thetaCoordinate;

    if (ri < 0) && (HankelType == 2)
      temporalRay.hankelType = 1; % change hankel type of ray
      ri= abs(ri);
      if (thi<0)
        thi = pi-abs(thi);
      else
        thi = abs(thi)-pi;
      end
      temporalRay.rCoordinate     = ri;
      temporalRay.thetaCoordinate = thi;
    end %conditions for cross origin

    
    % calculating Hankels
    HLr  = HankelLaguerre(r ,thi,LParametersZi,temporalRay.hankelType);
    HLth = HankelLaguerre(ri,th ,LParametersZi,temporalRay.hankelType);
    HLz  = HankelLaguerre(ri,thi,LParametersZ ,temporalRay.hankelType);     

    fr   = unwrap(angle(HLr.OpticalFieldLaguerre));
    fth  = unwrap(angle(HLth.OpticalFieldLaguerre));
    fz   = unwrap(angle(HLz.OpticalFieldLaguerre));

    % Calculating gradient
    [temporalRay] = gradientCylindrical(fr,fth,fz,LParametersZi.k,difr,temporalRay);

    % copying new coordinates of ray to object that includes all rays
    Rays = copyRay2ArrayRay(temporalRay,Rays,point_index);                                
  end

end