classdef GaussianParameters
    % GaussianParameters - Gaussian beam parameters
    % Compatible with GNU Octave and MATLAB
    %
    % This class serves two roles:
    %
    %   1. Snapshot at a fixed z (original API, fully preserved):
    %      params = GaussianParameters(z, w0, lambda)
    %      params.Waist, params.GouyPhase, params.Radius  -- scalar/vector at z
    %
    %   2. Dynamic evaluation at any z (new API, Phase 2):
    %      params = GaussianParameters(z, w0, lambda)  -- z can be 0 or any reference
    %      w   = params.waist(z2)       -- waist at arbitrary z2
    %      psi = params.gouyPhase(z2)   -- Gouy phase at arbitrary z2
    %      R   = params.radius(z2)      -- radius of curvature at arbitrary z2
    %      a   = params.amplitude(z2)   -- 1/waist at arbitrary z2
    %      al  = params.computeAlphaAtZ(z2)  -- complex beam parameter alpha at z2
    %
    % Both APIs are correct; propagators and beam methods should prefer the
    % dynamic API (role 2) so that parameters are not tied to a single z.

    properties
        zCoordinate
        InitialWaist
        Lambda          % Wavelength (m); alias Wavelength for backward compat.
        RayleighDistance
        k
        Waist
        GouyPhase
        Radius
        Amplitude
        DivergenceAngle
    end

    properties (Dependent)
        Wavelength      % Alias for Lambda (backward compatibility)
    end

    methods
        function obj = GaussianParameters(z, w0, lambda)
            if nargin > 0
                obj.zCoordinate = z;
                obj.InitialWaist = w0;
                obj.Lambda = lambda;

                obj.RayleighDistance = PhysicalConstants.rayleighDistance(w0, lambda);
                obj.k = PhysicalConstants.waveNumber(lambda);

                obj.Waist          = PhysicalConstants.waistAtZ(w0, z, lambda, obj.RayleighDistance);
                obj.GouyPhase      = PhysicalConstants.gouyPhase(z, obj.RayleighDistance);
                obj.Radius         = PhysicalConstants.radiusOfCurvature(z, obj.RayleighDistance);
                obj.Amplitude      = 1 ./ obj.Waist;
                obj.DivergenceAngle = atan(w0 ./ obj.RayleighDistance);
            end
        end

        % -----------------------------------------------------------------
        % Dynamic evaluation methods (Phase 2 API)
        % -----------------------------------------------------------------

        function w = waist(obj, z)
            % waist  Beam waist w(z) = w0 * sqrt(1 + (z/zR)^2).
            %
            % Unlike the stored property Waist (fixed at construction z),
            % this method evaluates at any scalar or vector z.
            w = PhysicalConstants.waistAtZ(obj.InitialWaist, z, obj.Lambda, obj.RayleighDistance);
        end

        function psi = gouyPhase(obj, z)
            % gouyPhase  Gouy phase psi(z) = arctan(z / zR).
            %
            % Returns the cumulative Gouy phase at axial position z.
            psi = PhysicalConstants.gouyPhase(z, obj.RayleighDistance);
        end

        function R = radius(obj, z)
            % radius  Radius of curvature R(z) = z*(1 + (zR/z)^2).
            %
            % Returns Inf at z = 0 (plane wavefront at the waist).
            R = PhysicalConstants.radiusOfCurvature(z, obj.RayleighDistance);
        end

        function a = amplitude(obj, z)
            % amplitude  On-axis amplitude = 1/w(z).
            %
            % Proportional to the peak field amplitude of the Gaussian mode.
            a = 1 ./ obj.waist(z);
        end

        function a = computeAlphaAtZ(obj, z)
            % computeAlphaAtZ  Complex beam parameter alpha at arbitrary z.
            %
            % alpha(z) = i*k / (2*q(z))  where  q(z) = z + i*zR
            %
            % Used by Elegant Hermite/Laguerre variants. This method
            % evaluates alpha at the given z (unlike computeAlpha which
            % uses the stored zCoordinate).
            q = z + 1i * obj.RayleighDistance;
            a = 1i * obj.k / (2 * q);
        end

        % -----------------------------------------------------------------
        % Original methods (unchanged)
        % -----------------------------------------------------------------

        function str = toString(obj)
            str = sprintf(...
                'GaussianParameters:\n  zCoordinate: %g\n  InitialWaist: %g\n  Lambda: %g\n  RayleighDistance: %g\n  k: %g\n  Waist: %g\n', ...
                obj.zCoordinate(1), obj.InitialWaist, obj.Lambda, obj.RayleighDistance, obj.k, obj.Waist(1));
        end

        function res = isEqual(obj, other)
            res = abs(obj.zCoordinate - other.zCoordinate) < 1e-12 && ...
                  abs(obj.InitialWaist - other.InitialWaist) < 1e-12 && ...
                  abs(obj.Lambda - other.Lambda) < 1e-12;
            res = all(res);
        end

        function val = get.Wavelength(obj)
            val = obj.Lambda;
        end

        function obj = set.Wavelength(obj, val)
            obj.Lambda = val;
        end

        function a = computeAlpha(obj)
            % computeAlpha  Complex beam parameter at the stored zCoordinate.
            %
            % alpha = i*k / (2*q(z))  where  q(z) = z + i*zR
            %
            % Prefer computeAlphaAtZ(z) for dynamic evaluation at any z.
            q = obj.zCoordinate + 1i * obj.RayleighDistance;
            a = 1i * obj.k / (2 * q);
        end
    end

    methods (Static)
        function w = getWaist(z, w0, zr)
            w = w0 .* sqrt(1 + (z./zr).^2);
        end

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
