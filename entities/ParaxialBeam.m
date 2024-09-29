% entities/ParaxialBeam.m
classdef ParaxialBeam < IBeam
    properties
        wavelength  % Wavelength of the beam (meters)
        waist       % Beam waist at the origin (meters)
        position    % Starting position (z = 0)
    end
    
    methods
        function obj = ParaxialBeam(wavelength, waist, position)
            obj.wavelength = wavelength;
            obj.waist = waist;
            obj.position = position;
        end
        
        function z_r = rayleigh_range(obj)
            % Compute the Rayleigh range
            z_r = pi * obj.waist^2 / obj.wavelength;
        end
        
        function w = beam_waist(obj, z)
            % Compute the beam waist at distance z
            z_r = obj.rayleigh_range();
            w = obj.waist * sqrt(1 + (z / z_r)^2);
        end
        
        function R = radius_of_curvature(obj, z)
            % Compute the radius of curvature of the wavefront at distance z
            z_r = obj.rayleigh_range();
            if z == 0
                R = Inf;  % At the waist, curvature is infinite
            else
                R = z * (1 + (z_r / z)^2);
            end
        end
        
        function phi = phase_shift(obj, z)
            % Compute the phase shift at distance z
            z_r = obj.rayleigh_range();
            phi = atan(z / z_r);
        end
        
        function I = intensity_profile(obj, z, x, y)
            % Compute the intensity profile at (x, y, z)
            w_z = obj.beam_waist(z);
            I0 = 2 / (pi * obj.waist^2);  % Peak intensity at the beam waist
            r2 = x.^2 + y.^2;  % Radial distance squared
            I = I0 * (obj.waist / w_z)^2 * exp(-2 * r2 / w_z^2);
        end
    end
end
