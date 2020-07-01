classdef ParaxialBeams
  properties
    RayleighDistance
    InitialWaist
    PropagationDistance
    Waist
    Radius
    Phase
    Amplitude
    Type
  end
  
  methods
    
    function beam = ParaxialBeams(RayleighDistance,InitialWaist,PropagationDistance)
      if nargin > 0 
        
        beam.InitialWaist        = InitialWaist;
        beam.RayleighDistance    = RayleighDistance;
        beam.PropagationDistance = PropagationDistance;
        beam.Waist               = waistPhysicalGaussianBeam(PropagationDistance,InitialWaist,RayleighDistance);
        beam.Radius              = radiusPhysicalGaussianBeam(PropagationDistance,RayleighDistance);
        beam.Phase               = phasePhysicalGaussianBeam(PropagationDistance,RayleighDistance);
      end 
      
    end
    
  end
  
end