% services/BeamService.m
classdef BeamService
    properties
        beam  % IBeam interface
    end
    
    methods
        function obj = BeamService(beam)
            obj.beam = beam;  % Dependency injection of any beam type implementing IBeam
        end
        
        function waist = calculateBeamWaist(obj, z)
            waist = obj.beam.beam_waist(z);
        end
        
        function intensity = calculateIntensity(obj, z, x, y)
            intensity = obj.beam.intensity_profile(z, x, y);
        end
        
        function phase = calculatePhaseShift(obj, z)
            phase = obj.beam.phase_shift(z);
        end
    end
end
