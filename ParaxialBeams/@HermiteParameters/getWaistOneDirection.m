function waistH = getWaistOnedirection(zDistance,InitialWaist,RayleighDistance,l)

  waistH = sqrt(2*l+1)...
         .*(InitialWaist).*sqrt((zDistance./RayleighDistance).^2+1); % factor of Gaussian waist

end