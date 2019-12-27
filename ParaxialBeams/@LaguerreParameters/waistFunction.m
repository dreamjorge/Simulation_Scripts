function LWaistZ = waistFunction(distancePropagation,InitialWaist,RayleighDistance,nu,mu)
%This function estimates waist of Laguerre Gaussian Beam
% receives distancePropagation,initialWaist,RayleighDistance as inputs and
% gives waist of Gaussian Beam
%
%Example: LWaistZ = waistFunction(distancePropagation,InitialWaist,RayleighDistance,nu,mu)

  LWaistZ = sqrt(2*nu+mu+1)...
          .*(InitialWaist).*sqrt((distancePropagation./RayleighDistance).^2+1); % factor of Gaussian waist

end