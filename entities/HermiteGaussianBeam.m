% entities/HermiteGaussianBeam.m
classdef HermiteGaussianBeam < ParaxialBeam
    properties
        m  % Mode index in x
        n  % Mode index in y
    end
    
    methods
        function obj = HermiteGaussianBeam(wavelength, waist, position, m, n)
            obj@ParaxialBeam(wavelength, waist, position);  % Call the base constructor
            obj.m = m;
            obj.n = n;
        end
        
        function I = intensity_profile(obj, z, x, y)
            % Intensity profile for Hermite-Gaussian beams
            w_z = obj.beam_waist(z);
            Hm = hermiteH(obj.m, sqrt(2) * x / w_z);
            Hn = hermiteH(obj.n, sqrt(2) * y / w_z);
            I0 = 2 / (pi * obj.waist^2);
            I = I0 * (obj.waist / w_z)^2 * Hm.^2 .* Hn.^2 .* exp(-2 * (x.^2 + y.^2) / w_z^2);
        end
        
        function phi = phase_shift(obj, z)
            % Hermite-Gaussian phase shift adds an extra mode term
            phi = (obj.m + obj.n + 1) * atan(z / obj.rayleigh_range());
        end
    end
end
