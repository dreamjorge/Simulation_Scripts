function [Rays] = getPropagateCartesianRays(Rays,...
                                            x,y,...
                                            difr,...
                                            HParametersZi,...
                                            HParametersZ,...
                                            HankelType) 

%% Obtain New coordinates of ray using a step of ray in z-direction,
%equation or ray r(z) = mrz*dz+r(z-1) using previous value
%debugMode = 'Off';                                                         
  % each component of diferential dr 
  tempdr    = num2cell(difr);
  [~,~,dz]  = deal(tempdr{:});
  k         = HParametersZi.k;
  TotalRays = numel(Rays.xCoordinate);
  
  for ray_index = 1 : TotalRays

    temporalRay = OpticalRay(); 

    % copy data a temporal ray 
    temporalRay = copyArrayRay2Ray(Rays,temporalRay,ray_index);
    
    % new coordinates of ray
    temporalRay.xCoordinate = temporalRay.xCoordinate + (1./temporalRay.zxSlope)*dz;
    temporalRay.yCoordinate = temporalRay.yCoordinate + (1./temporalRay.zySlope)*dz;
    temporalRay.zCoordinate = temporalRay.zCoordinate + dz;
  
    xi   = temporalRay.xCoordinate;
    yi   = temporalRay.yCoordinate;
    
    % calculating Hankels
    HHx  = HankeleHermite(x ,yi,HParametersZi,HankelType);
    HHy  = HankeleHermite(xi,y ,HParametersZi,HankelType);
    HHz  = HankeleHermite(xi,yi,HParametersZ ,HankelType);     

    fx   = unwrap(angle(HHx.OpticalField));
    fy   = unwrap(angle(HHy.OpticalField));
    fz   = unwrap(angle(HHz.OpticalField));

    % Calculating gradient
    [temporalRay] = getCartesianGradient(fx,fy,fz,k,difr,temporalRay);
    
    % copying new coordinates of ray to object that includes all rays
    Rays = copyRay2ArrayRay(temporalRay,Rays,ray_index);                                
  end

end