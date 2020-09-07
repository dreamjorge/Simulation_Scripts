function [Rays] = getPropagateCylindricalRays(Rays,...
                                              TotalRays,...
                                              r,th,...
                                              difr,...
                                              LParametersZi,...
                                              LParametersZ...
                                              ) 
%% Obtain New coordinates of ray using a step of ray in z-direction,
%equation or ray r(z) = mrz*dz+r(z-1) using previous value
                                                         
  % each component of diferential dr 
  tempdr     = num2cell(difr);
  [~,~,dz] = deal(tempdr{:});
  
  temporalRay = CylindricalRay(); 

  for point_index = 1 : TotalRays
    
     % copy data a temporal ray 
     
     temporalRay = copyArrayRay2Ray(Rays,temporalRay,point_index);

     % new coordinates of ray

    temporalRay.rCoordinate     = temporalRay.rCoordinate     + (1./temporalRay.zrSlope )*dz;
    temporalRay.thetaCoordinate = temporalRay.thetaCoordinate + (1./temporalRay.zthSlope)*dz;
    temporalRay.zCoordinate     = temporalRay.zCoordinate + dz;
    % obtain new (x,y) coordinates 
    [Rays.xCoordinate,Rays.yCoordinate] = pol2cart(Rays.thetaCoordinate,Rays.rCoordinate);

    % taking this new coordinates of point to calculate hankels
    ri  = temporalRay.rCoordinate;
    thi = temporalRay.thetaCoordinate;

    if (ri < 0) && (temporalRay.hankelType == 2)
      temporalRay.hankelType = 1;
    end 

    temporalRay.rCoordinate = abs(ri);
    
    % calculating Hankels
    HLr  = HankelLaguerre(r ,thi,LParametersZi,temporalRay.hankelType);
    HLth = HankelLaguerre(ri,th ,LParametersZi,temporalRay.hankelType);
    HLz  = HankelLaguerre(ri,thi,LParametersZ ,temporalRay.hankelType);     

    fr   = unwrap(angle(HLr.OpticalFieldLaguerre));
    fth  = unwrap(angle(HLth.OpticalFieldLaguerre));
    fz   = unwrap(angle(HLz.OpticalFieldLaguerre));

   % Calculating gradient
    [temporalRay] = getGradientCylindrical(fr,fth,fz,LParametersZi.k,difr,temporalRay);
    
    if (temporalRay.hankelType == 1)
      temporalRay.zrSlope = abs(temporalRay.zrSlope);
    elseif (temporalRay.hankelType == 2)
      temporalRay.zrSlope = -abs(temporalRay.zrSlope);
    end 
    % copying new coordinates of ray to object that includes all rays
    Rays = copyRay2ArrayRay(temporalRay,Rays,point_index);                                
  end    
                    


end