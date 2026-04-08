classdef GaussianParameters
    % GaussianParameters - Gaussian beam parameters
    % Compatible with GNU Octave and MATLAB
    
    properties
        zCoordinate
        InitialWaist
        Wavelength
        RayleighDistance
        k
        Waist
        GouyPhase
        Radius
        Amplitude
        DivergenceAngle
    end
    
    methods
        function obj = GaussianParameters(z, w0, lambda)
            if nargin > 0
                obj.zCoordinate = z;
                obj.InitialWaist = w0;
                obj.Wavelength = lambda;
                
                % Use PhysicalConstants for calculations
                obj.RayleighDistance = PhysicalConstants.rayleighDistance(w0, lambda);
                obj.k = PhysicalConstants.waveNumber(lambda);
                
                % Element-wise calculations for vector z
                obj.Waist = PhysicalConstants.waistAtZ(w0, z, lambda, obj.RayleighDistance);
                obj.GouyPhase = PhysicalConstants.gouyPhase(z, obj.RayleighDistance);
                obj.Radius = PhysicalConstants.radiusOfCurvature(z, obj.RayleighDistance);
                obj.Amplitude = 1 ./ obj.Waist;
                obj.DivergenceAngle = atan(w0 ./ obj.RayleighDistance);
            end
        end
        
        function str = toString(obj)
            str = sprintf(...
                'GaussianParameters:\n  zCoordinate: %g\n  InitialWaist: %g\n  Wavelength: %g\n  RayleighDistance: %g\n  k: %g\n  Waist: %g\n', ...
                obj.zCoordinate(1), obj.InitialWaist, obj.Wavelength, obj.RayleighDistance, obj.k, obj.Waist(1));
        end
        
        function res = isEqual(obj, other)
            res = abs(obj.zCoordinate - other.zCoordinate) < 1e-12 && ...
                  abs(obj.InitialWaist - other.InitialWaist) < 1e-12 && ...
                  abs(obj.Wavelength - other.Wavelength) < 1e-12;
            res = all(res);
        end
    end
    
    methods
        function a = computeAlpha(obj)
            % Complex beam parameter: alpha = i*k / (2*q(z))
            % where q(z) = z + i*zR (complex beam parameter)
            % Used by ElegantHermiteParameters and ElegantLaguerreParameters.
            q = obj.zCoordinate + 1i * obj.RayleighDistance;
            a = 1i * obj.k / (2 * q);
        end
    end

    methods (Static)
        function w = getWaist(z, w0, zr)
            w = w0 .* sqrt(1 + (z./zr).^2);
        end

        % These delegate to PhysicalConstants for backward compatibility.
        function zr = rayleighDistance(w0, lambda)
            zr = PhysicalConstants.rayleighDistance(w0, lambda);
        end

        function phase = getPhase(z, zr)
            phase = PhysicalConstants.gouyPhase(z, zr);
        end

        function R = getRadius(z, zr)
            R = PhysicalConstants.radiusOfCurvature(z, zr);
        end
    end
end
