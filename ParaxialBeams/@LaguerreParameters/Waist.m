function Lwaist = Waist(zDistance,InitialWaist,RayleighDistance,nu,mu)
%This function estimates waist of Laguerre Gaussian Beam
% receives distancePropagation,initialWaist,RayleighDistance, l, p  as inputs and
% gives waist of Gaussian Beam
  Lwaist = sqrt(2)*sqrt(2*nu+mu+1)...
         .*(InitialWaist).*sqrt((zDistance./RayleighDistance).^2+1); % factor of Gaussian waist

end