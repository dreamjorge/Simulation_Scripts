function [] = plotLaguerreParameters(LGBParameters)
%% Function plots parameters of Gaussian Beam (Waist,DivergenceAngle)
% Input:
%  -GB Parameters as GaussianBeamParameters

p1            = plot(LGBParameters.zCoordinate,LGBParameters.Waist,'Color','red');
hold on
p2            = plot(LGBParameters.zCoordinate,-LGBParameters.Waist,'Color','red');
p3            = plot(LGBParameters.zCoordinate, LGBParameters.zCoordinate*tan(LGBParameters.DivergenceAngle),'Color','blue');
p4            = plot(LGBParameters.zCoordinate,-LGBParameters.zCoordinate*tan(LGBParameters.DivergenceAngle),'Color','blue');
p5            = plot(LGBParameters.zCoordinate,-LGBParameters.LaguerreWaist,'Color','green');
p6            = plot(LGBParameters.zCoordinate, LGBParameters.LaguerreWaist,'Color','green');
xlabel('Distance of Propagation [microns]')
ylabel('[microns]')
title('Parameters of Laguerre Gaussian Beam')
legend([p1,p5,p3],{'Waist of Gaussian Beam'...
                  ,'Waist of Laguerre Beam'...
                  ,['Angle of Divergence = ',num2str(rad2deg(LGBParameters.DivergenceAngle)),'°']})
hold off
end