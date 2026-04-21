classdef BeamComputation
    % BeamComputation - Stateless computation utilities for beam parameters
    % Compatible with GNU Octave and MATLAB
    %
    % Purpose:
    %   Provides pure, stateless static methods for Gaussian beam formulas.
    %   No object state, no side effects. Same inputs -> same outputs.
    %
    % Usage:
    %   zr  = BeamComputation.rayleighDistance(w0, lambda);
    %   k   = BeamComputation.waveNumber(lambda);
    %   w   = BeamComputation.waist(w0, z, lambda, zr);
    %   psi = BeamComputation.gouyPhase(z, zr);
    %   R   = BeamComputation.radiusOfCurvature(z, zr);
    %   q   = BeamComputation.complexBeamParameter(z, zr, k);
    %
    % Notes:
    %   Formulas match established optics references (Kogelnik & Li, Siegman).
    %   All methods are fully vectorized.

    methods (Static)

        function zr = rayleighDistance(w0, lambda)
            % rayleighDistance  Rayleigh range zR = pi * w0^2 / lambda
            %
            % Input:
            %   w0     beam waist at z=0 (m)
            %   lambda wavelength (m)
            %
            % Output:
            %   zr     Rayleigh range (m)

            zr = pi * w0.^2 ./ lambda;
        end


        function k = waveNumber(lambda)
            % waveNumber  Wave number k = 2*pi/lambda
            %
            % Input:
            %   lambda wavelength (m)
            %
            % Output:
            %   k      wave number (rad/m)

            k = 2 * pi ./ lambda;
        end


        function w = waist(w0, z, lambda, zr)
            % waist  Beam waist w(z) = w0 * sqrt(1 + (z/zR)^2)
            %
            % Input:
            %   w0     beam waist at z=0 (m)
            %   z      axial position (m), scalar or vector
            %   lambda wavelength (m)
            %   zr     Rayleigh range (m)
            %
            % Output:
            %   w      beam waist at z (m)

            if nargin < 4
                zr = BeamComputation.rayleighDistance(w0, lambda);
            end
            w = w0 .* sqrt(1 + (z ./ zr).^2);
        end


        function psi = gouyPhase(z, zr)
            % gouyPhase  Gouy phase psi(z) = arctan(z / zR)
            %
            % Input:
            %   z      axial position (m), scalar or vector
            %   zr     Rayleigh range (m)
            %
            % Output:
            %   psi    Gouy phase (rad)

            psi = atan(z ./ zr);
        end


        function R = radiusOfCurvature(z, zr)
            % radiusOfCurvature  Radius of curvature R(z) = z*(1 + (zR/z)^2)
            %
            % At z=0 the wavefront is planar: R -> Inf
            %
            % Input:
            %   z      axial position (m), scalar or vector
            %   zr     Rayleigh range (m)
            %
            % Output:
            %   R      radius of curvature (m), Inf at z=0

            R = z .* (1 + (zr ./ z).^2);
            R(z == 0) = Inf;
        end


        function q = complexBeamParameter(z, zr, k)
            % complexBeamParameter  Complex beam parameter q(z) = z + i*zR
            %
            % The complex beam parameter is central to "elegant" beam families.
            % Its reciprocal 1/q(z) appears in the argument of Hermite/Laguerre
            % polynomials for elegant variants.
            %
            % Input:
            %   z      axial position (m)
            %   zr     Rayleigh range (m)
            %   k      wave number (rad/m)
            %
            % Output:
            %   q      complex beam parameter q(z) = z + i*zR

            q = z + 1i * zr;
        end


        function alpha = complexAlpha(z, zr, k)
            % complexAlpha  Complex alpha = i*k/(2*q(z))
            %
            % Used by ElegantHermiteBeam and ElegantLaguerreBeam as the
            % scale factor in the complex argument to Hermite/Laguerre polynomials.
            %
            % Input:
            %   z      axial position (m)
            %   zr     Rayleigh range (m)
            %   k      wave number (rad/m)
            %
            % Output:
            %   alpha  complex beam parameter alpha (1/m)

            q = z + 1i * zr;
            alpha = 1i * k ./ (2 .* q);
        end

    end
end
