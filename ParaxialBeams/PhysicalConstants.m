classdef PhysicalConstants
    % PhysicalConstants - Physical constants and utility functions
    % Compatible with GNU Octave and MATLAB
    
    properties (Constant)
        speed_of_light = 299792458;
        planck = 6.62607015e-34;
        planck_reduced = 1.054571817e-34;
        vacuum_permittivity = 8.8541878128e-12;
        vacuum_permeability = 1.25663706212e-6;
        impedance_vacuum = 376.730313668;
    end
    
    methods (Static)
        function k = waveNumber(lambda)
            k = 2*pi ./ lambda;
        end
        
        function zr = rayleighDistance(w0, lambda)
            zr = pi * w0.^2 ./ lambda;
        end
        
        function w = waistAtZ(w0, z, lambda, zr)
            if nargin < 4
                zr = PhysicalConstants.rayleighDistance(w0, lambda);
            end
            w = w0 .* sqrt(1 + (z ./ zr).^2);
        end
        
        function R = radiusOfCurvature(z, zr)
            % At z=0 the wavefront is flat: R -> Inf (avoids 0*Inf = NaN)
            R = z .* (1 + (zr ./ z).^2);
            R(z == 0) = Inf;
        end
        
        function gouy = gouyPhase(z, zr)
            gouy = atan(z ./ zr);
        end
    end
end
