function [Rays] = getPropagateRay(Rays,...
                                  TotalRays,...
                                  r,th,...
                                  difr,...
                                  LParametersZi,...
                                  LParametersZ,...
                                  HankelType) 
                                
if HankelType == 2
  tempHankel = 2;
else
  tempHankel = 1;
end
%% Obtain New coordinates of ray using a step of ray in z-direction,
%equation or ray r(z) = mrz*dz+r(z-1) using previous value
debugMode = 'Off';                                                         
  % each component of diferential dr 
  tempdr     = num2cell(difr);
  [dr,~,dz] = deal(tempdr{:});
 
  for point_index = 1 : TotalRays

    temporalRay = CylindricalRay(); 

    % copy data a temporal ray 
    temporalRay = copyArrayRay2Ray(Rays,temporalRay,point_index);
    % new coordinates of ray
    rtemp  = temporalRay.rCoordinate ;
    thtemp = temporalRay.thetaCoordinate;
    temporalRay.rCoordinate     = temporalRay.rCoordinate     + (1./temporalRay.zrSlope )*dz;
    temporalRay.thetaCoordinate = temporalRay.thetaCoordinate + (1./temporalRay.zthSlope)*dz;
    temporalRay.zCoordinate     = temporalRay.zCoordinate + dz;
    % obtain new (x,y) coordinates 
    [Rays.xCoordinate,Rays.yCoordinate] = pol2cart(Rays.thetaCoordinate,Rays.rCoordinate);

    % taking this new coordinates of point to calculate hankels
    ri  = temporalRay.rCoordinate;
    thi = temporalRay.thetaCoordinate;

    if (ri < 0) && (tempHankel == 2)
      tempHankel = 1;
    end 

    temporalRay.rCoordinate = abs(ri);
    
    % calculating Hankels
    HLr  = HankelLaguerre(r ,thi,LParametersZi,tempHankel);
    HLth = HankelLaguerre(ri,th ,LParametersZi,tempHankel);
    HLz  = HankelLaguerre(ri,thi,LParametersZ ,tempHankel);     

    fr   = unwrap(angle(HLr.OpticalFieldLaguerre));
    fth  = unwrap(angle(HLth.OpticalFieldLaguerre));
    fz   = unwrap(angle(HLz.OpticalFieldLaguerre));

if strcmp(debugMode,'On')
    if (HankelType == 1)
    figure(100)
    subplot(3,2,point_index)
    plot(r,fr)
    vline(ri)
    title('Hr')
    xlabel('r')
    subplot(3,2,point_index+2)
    plot(th,fth)
    vline(thi)
    title('Hth')
    xlabel('th')
    subplot(3,2,point_index+4)
    plot(LParametersZ.zCoordinate,fz)
    vline(LParametersZi.zCoordinate)
    title('Hz')
    xlabel('z')
    suptitle('Hankel 1')
    elseif(HankelType == 2)
    gr = gradient(fr)/dr;
    N  = size(gr,2);
    slr=gr(N/2+1+floor(ri/dr));
    if (point_index == 2)
      if slr >0
        print(num2str(LParametersZi.zCoordinate))
      end
    end
    figure(101)
    suptitle('Hankel 2')
    subplot(3,2,point_index)
    plot(r,fr)
    vline(ri)
    title('Hr')
    xlabel(['r'])
    subplot(3,2,point_index+2)
    plot(th,fth)
    vline(thi)
    title('Hth')
    xlabel('th')
    subplot(3,2,point_index+4)
    plot(LParametersZ.zCoordinate,fz)
    vline(LParametersZi.zCoordinate)
    title('Hz')
    xlabel('z')

    %suptitle(['phases at point (',num2str(ri),',',num2str(thi),',',num2str(temporalRay.zCoordinate),')'])
    end
end
    % Calculating gradient
    [temporalRay] = gradientCylindrical(fr,fth,fz,LParametersZi.k,difr,temporalRay);
    
    if (tempHankel == 1)
      temporalRay.zrSlope = abs(temporalRay.zrSlope);
    elseif (tempHankel == 2)
      temporalRay.zrSlope = -abs(temporalRay.zrSlope);
    end 
    % copying new coordinates of ray to object that includes all rays
    Rays = copyRay2ArrayRay(temporalRay,Rays,point_index);                                
  end

end