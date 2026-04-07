classdef PhysicalConstants
%PHYSICALCONSTANTS Physical constants for optical simulations
%   This class provides fundamental and derived constants for beam
%   propagation simulations.
%
%   Example:
%       c = PhysicalConstants;
%       k = c.hbar * c.k; % Planck-reduced wave number
    
    properties (Constant)
        %% Fundamental constants (SI units)
        c           speed_of_light         % Speed of light in vacuum (m/s)
        h           planck                  % Planck constant (J·s)
        hbar        planck_reduced          % Reduced Planck constant (J·s)
        e           elementary_charge       % Elementary charge (C)
        eps0        vacuum_permittivity     % Vacuum permittivity (F/m)
        mu0         vacuum_permeability     % Vacuum permeability (H/m)
        
        %% Derived optical constants
        Z0          impedance_vacuum        % Vacuum impedance (Ω)
    end
    
    methods (Static)
        function self = PhysicalConstants()
            % Constructor - constants only, no instance needed
            if nargin > 0
                error('PhysicalConstants is a static class');
            end
        end
    end
    
    %% Get methods for fundamental constants
    methods (Static)
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
            val = 376.730313668; % sqrt(mu0/eps0)
        end
    end
    
    %% Utility methods
    methods (Static)
        function k = waveNumber(wavelength)
            %WAVENUMBER Calculate wave number k = 2*pi/lambda
            %   k = waveNumber(lambda) returns k in 1/meters
            %
            %   Example:
            %       k = waveNumber(632.8e-9); % HeNe laser
            
            k = 2*pi ./ wavelength;
        end
        
        function zr = rayleighDistance(w0, wavelength)
            %RAYLEIGH Distance Calculate Rayleigh distance
            %   zr = rayleighDistance(w0, lambda) returns Rayleigh range
            %   zr = pi*w0^2/lambda
            %
            %   Example:
            %       zr = rayleighDistance(100e-6, 632.8e-9);
            
            zr = pi * w0.^2 ./ wavelength;
        end
        
        function w = waistAtZ(w0, z, wavelength, zr)
            %WAISTATZ Calculate beam waist at propagation distance
            %   w = waistAtZ(w0, z, lambda, zr) returns beam radius at z
            %   w = w0*sqrt(1 + (z/zr)^2)
            %
            %   Example:
            %       w = waistAtZ(100e-6, 0.1, 632.8e-9, 0.05);
            
            if nargin < 4
                zr = PhysicalConstants.rayleighDistance(w0, wavelength);
            end
            w = w0 .* sqrt(1 + (z ./ zr).^2);
        end
        
        function R = radiusOfCurvature(z, wavelength, zr)
            %RADIUSOFCURVATURE Calculate radius of curvature
            %   R = radiusOfCurvature(z, lambda, zr) returns R = z*(1 + (zr/z)^2)
            %
            %   Example:
            %       R = radiusOfCurvature(0.05, 632.8e-9, 0.05);
            
            R = z .* (1 + (zr ./ z).^2);
        end
        
        function gouy = gouyPhase(z, zr)
            %GOUYPHASE Calculate Gouy phase
            %   gouy = gouyPhase(z, zr) returns arctan(z/zr)
            %
            %   Example:
            %       gouy = gouyPhase(0.05, 0.1);
            
            gouy = atan(z ./ zr);
        end
    end
end
