% entities/IBeam.m
classdef (Abstract) IBeam
    methods (Abstract)
        beam_waist(obj, z)           % Calculate the beam waist at position z
        radius_of_curvature(obj, z)  % Calculate the wavefront's radius of curvature
        phase_shift(obj, z)          % Compute the phase shift at position z
        intensity_profile(obj, z, x, y)  % Get the intensity at position (x, y)
    end
end
