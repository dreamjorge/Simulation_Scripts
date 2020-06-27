function Lwaist = getWaist(zDistance,InitialWaist,RayleighDistance,nu,mu)

  Lwaist = sqrt(2)*sqrt(2*nu+mu+1)...
         .*(InitialWaist).*sqrt((zDistance./RayleighDistance).^2+1); % factor of Gaussian waist

end