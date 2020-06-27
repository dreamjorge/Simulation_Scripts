function waistL = getWaist(zDistance,InitialWaist,RayleighDistance,l)

  waistL = sqrt(2*l+1)...
         .*(InitialWaist).*sqrt((zDistance./RayleighDistance).^2+1); % factor of Gaussian waist

end