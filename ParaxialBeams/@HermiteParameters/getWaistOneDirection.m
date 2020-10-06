function waistH = getWaistOneDirection(zDistance,InitialWaist,RayleighDistance,n)
%% Function for Obtain of waist of Hermite Guassian in 1-D
% Recives next elements
% - zDistance         Distance in z for estimate waist.
% - InitialWaist      Waist of Gaussian Beam at 0 distance.
% - RayleighDistance  Distance of Rayleigh of Guassian Beam.
% - n                 Number of Poylinolial Hermite.

  waistH = sqrt(2)*sqrt(2*n+1)...
         .*(InitialWaist).*sqrt((zDistance./RayleighDistance).^2+1); % factor of Gaussian waist

end