function radiusZ = radiusFunction(distancePropagation,RayleighDistance)
%This function estimates radius of Gaussian Beam
% receives distancePropagation,initialWaist,RayleighDistance as inputs and
% gives phase of Gaussian Beam
%
%Example: radiusZ = radiusFunction(distancePropagation,RayleighDistance)
  radiusZ = (distancePropagation).*(1+(RayleighDistance./distancePropagation).^2);

  if (find(isnan(radiusZ))~= 0)
      radiusZ(isnan(radiusZ))= inf;
  end
end