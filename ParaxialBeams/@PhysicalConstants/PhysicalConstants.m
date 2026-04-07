classdef PhysicalConstants
%PHYSICALCONSTANTS Physical constants for optical simulations
%   This class provides fundamental and derived constants for beam
%   propagation simulations.
    
    properties (Constant)
        c           speed_of_light
        h           planck
        hbar        planck_reduced
        e           elementary_charge
        eps0        vacuum_permittivity
        mu0         vacuum_permeability
        Z0          impedance_vacuum
    end
    
    methods (Static)
        function self = PhysicalConstants()
            if nargin > 0
                error('PhysicalConstants is a static class');
            end
        end
        
        function val = speed_of_light()
            val = 299792458;
        end
        
        function val = planck()
            val = 6.62607015e-34;
        end
        
        function val = planck_reduced()
            val = 1.054571817e-34;
        end
        
        function val = elementary_charge()
            val = 1.602176634e-19;
        end
        
        function val = vacuum_permittivity()
            val = 8.8541878128e-12;
        end
        
        function val = vacuum_permeability()
            val = 1.25663706212e-6;
        end
        
        function val = impedance_vacuum()
            val = 376.730313668;
        end
        
        function k = waveNumber(wavelength)
            k = 2*pi ./ wavelength;
        end
        
        function zr = rayleighDistance(w0, wavelength)
            zr = pi * w0.^2 ./ wavelength;
        end
        
        function w = waistAtZ(w0, z, wavelength, zr)
            if nargin < 4
                zr = PhysicalConstants.rayleighDistance(w0, wavelength);
            end
            w = w0 .* sqrt(1 + (z ./ zr).^2);
        end
        
        function R = radiusOfCurvature(z, zr)
            R = z .* (1 + (zr ./ z).^2);
        end
        
        function gouy = gouyPhase(z, zr)
            gouy = atan(z ./ zr);
        end
    end
end
