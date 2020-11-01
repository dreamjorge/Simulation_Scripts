function PhiZ = getPhase(distancePropagation,RayleighDistance)
%This function estimates phase of Gaussian Beam
% receives distancePropagation,RayleighDistance as inputs and
% gives phase of Gaussian Beam
%
%Example: PhiZ = getPhase(distancePropagation,RayleighDistance)
  PhiZ = atan(distancePropagation/RayleighDistance);

end