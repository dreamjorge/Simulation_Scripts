function [Rays] = getPropagateCylindricalRays(Rays,...
                                              TotalRays,...
                                              r,th,...
                                              difr,...
                                              LParametersZi,...
                                              LParametersZ,...
                                              HankelType) 
%% Obtain New coordinates of ray using a step of ray in z-direction,
%equation or ray r(z) = mrz*dz+r(z-1) using previous value
                                                         
  % each component of diferential dr 
  tempdr     = num2cell(difr);
  [~,~,dz]   = deal(tempdr{:});
  k          = LParametersZi.k;
  
  temporalRay = CylindricalRay(); 

  for point_index = 1 : TotalRays
    
    % copy data a temporal ray 
    temporalRay = copyArrayRay2Ray(Rays,temporalRay,point_index);
    % slopes
    zrSlope                     = temporalRay.zrSlope;
    zthSlope                    = temporalRay.zthSlope;
    % propagate ray 
    temporalRay.rCoordinate     = temporalRay.rCoordinate     + (1./zrSlope )*dz;
    temporalRay.thetaCoordinate = temporalRay.thetaCoordinate + (1./zthSlope)*dz;
    temporalRay.zCoordinate     = temporalRay.zCoordinate     + dz;
   
    % taking this new coordinates of point to calculate hankels
    ri  = temporalRay.rCoordinate;
    thi = temporalRay.thetaCoordinate;

    
    % Condition for cross origin (Change hankel type) and keep r positive
    if (ri < 0) && (HankelType == 2)
      temporalRay.hankelType = 1;
      ri                     = abs(ri);
    end 
    
    temporalRay.rCoordinate     = ri;
    temporalRay.thetaCoordinate = thi;
    
    %if we have change of hankel we need change sign of x,y coordinates
    if((HankelType == 2) && (temporalRay.hankelType == 1))
      [xi,yi] = pol2cart(thi,ri);
      temporalRay.xCoordinate = - xi;
      temporalRay.yCoordinate = - yi;
      
    else
      [xi,yi] = pol2cart(thi,ri);
      temporalRay.xCoordinate = xi;
      temporalRay.yCoordinate = yi;
      
    end

    % calculating Hankels in point (r,th,z) of ray
    HLr  = HankelLaguerre(r ,thi,LParametersZi,temporalRay.hankelType);
    HLth = HankelLaguerre(ri,th ,LParametersZi,temporalRay.hankelType);
    HLz  = HankelLaguerre(ri,thi,LParametersZ ,temporalRay.hankelType);     

    fr   = unwrap(angle( HLr.OpticalFieldLaguerre));
    fth  = unwrap(angle(HLth.OpticalFieldLaguerre));
    fz   = unwrap(angle( HLz.OpticalFieldLaguerre));

     % Calculating gradient (avoid introduce r == 0 coordinate)
    [temporalRay] = getCylindricalGradient(fr,fth,fz,k,difr,temporalRay);
    
    % copying new coordinates of ray to object that includes all rays
    Rays = copyRay2ArrayRay(temporalRay,Rays,point_index);                                
  end    
                    


end