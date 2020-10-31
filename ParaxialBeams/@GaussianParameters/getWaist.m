function waistZ = getWaist(distancePropagation,initialWaist,RayleighDistance)
%This function estimates waist of Gaussian Beam
% receives distancePropagation,initialWaist,RayleighDistance as inputs and
% gives waist of Gaussian Beam
%
%Example: waistZ = getWaist(distancePropagation,initialWaist,RayleighDistance)
  waistZ = (initialWaist)*sqrt( (distancePropagation./RayleighDistance).^2 + 1);

end
