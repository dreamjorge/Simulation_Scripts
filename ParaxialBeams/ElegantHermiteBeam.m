classdef ElegantHermiteBeam < ParaxialBeam
    % ElegantHermiteBeam - Elegant Hermite-Gaussian beam implementation
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: CARTESIAN (x, y) - same as HermiteBeam.
    % The unified API opticalField(X, Y, z) is natively Cartesian.
    %
    % Mathematical differences from standard HermiteBeam:
    %
    %   Standard HG: H_n(sqrt(2)*x/w(z))    -- real argument, waist scaling
    %   Elegant HG:  H_n(sqrt(alpha(z))*x)   -- complex argument, alpha scaling
    %
    % where alpha(z) = i*k / (2*q(z)), q(z) = z + i*z_R (complex beam parameter).
    %
    % Full field definition (elegant Hermite-Gaussian, EHG_{nm}):
    %
    %   u_{nm}(x,y,z) = H_n(sqrt(alpha)*x) * H_m(sqrt(alpha)*y) *
    %                   u_0(r,z) * exp(i*(n+m)*psi(z))
    %
    % Physical consequence: because alpha is complex, the Hermite polynomials
    % evaluated at complex arguments produce amplitude AND phase modulation
    % simultaneously, which is the hallmark of the "elegant" variant.
    %
    % Reference: Siegman, "Lasers", University Science Books (1986), Ch. 17;
    %            Siegman, J. Opt. Soc. Am. A 13, 952 (1996).

    properties
        Parameters      % ElegantHermiteParameters object
        OpticalField    % 2D array representing the field (at constructor z)
        x               % x coordinate matrix (stored for reference)
        y               % y coordinate matrix (stored for reference)
    end

    methods
        function obj = ElegantHermiteBeam(x, y, params)
            % Constructor
            % x, y:   Cartesian coordinate matrices
            % params: ElegantHermiteParameters object

            if nargin > 0
                obj = obj@ParaxialBeam(params.Lambda);
                obj.Parameters   = params;
                obj.x            = x;
                obj.y            = y;
                obj.OpticalField = obj.computeField(x, y, params);
            else
                obj = obj@ParaxialBeam();
            end
        end

        % -----------------------------------------------------------------
        % ParaxialBeam interface
        % -----------------------------------------------------------------

        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex optical field on Cartesian grid (X,Y) at depth z.
            %
            % Reconstructs ElegantHermiteParameters at z (recomputes alpha(z))
            % before evaluating the field.
            w0       = obj.Parameters.InitialWaist;
            lambda   = obj.Lambda;
            n        = obj.Parameters.n;
            m        = obj.Parameters.m;
            params_z = ElegantHermiteParameters(z, w0, lambda, n, m);
            field    = obj.computeField(X, Y, params_z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            w0     = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            params = GaussianParameters(z, w0, lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'elegant_hermite_2_0'.
            name = sprintf('elegant_hermite_%d_%d', obj.Parameters.n, obj.Parameters.m);
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, x, y, params)
            % Compute EHG field from Cartesian grids and parameters.
            %
            % Key: sqrt(alpha) is complex, so H_n(sqrt(alpha)*x) is evaluated
            % with a complex argument, encoding both amplitude and phase.
            sqrt_alpha = sqrt(params.alpha);
            Hn = PolynomialUtils.hermitePoly(params.n, sqrt_alpha .* x);
            Hm = PolynomialUtils.hermitePoly(params.m, sqrt_alpha .* y);

            R  = sqrt(x.^2 + y.^2);
            GB = GaussianBeam(R, params);

            % PhiPhase = (n+m)*psi(z) — modal Gouy phase shift
            field = Hn .* Hm .* exp(1i * params.PhiPhase) .* GB.OpticalField;
        end
    end
end
