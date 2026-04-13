classdef HermiteBeam < ParaxialBeam
    % HermiteBeam - Scalar optical field for a Hermite-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: CARTESIAN (x, y)
    % Note: LaguerreBeam uses polar (r, theta). Both share the unified
    % opticalField(X, Y, z) API — coordinate conversion happens internally.
    %
    % Mathematical definition (standard Hermite-Gaussian, HG_{nm}):
    %
    %   u_{nm}(x,y,z) = H_n(sqrt(2)*x/w) * H_m(sqrt(2)*y/w) *
    %                   u_0(r,z) * exp(i*(n+m)*psi(z))
    %
    % where:
    %   H_n, H_m  - Hermite polynomials of order n, m
    %   w = w(z)  - beam waist at z (from GaussianParameters)
    %   u_0(r,z)  - fundamental Gaussian carrier field
    %   psi(z)    - Gouy phase = arctan(z / z_R)
    %
    % The extra Gouy phase factor exp(i*(n+m)*psi) is handled via
    % HermiteParameters.PhiPhase = (n+m)*psi(z).
    %
    % Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 3.

    properties
        Parameters      % HermiteParameters object
        OpticalField    % 2D array representing the field (at constructor z)
        x               % x coordinate matrix (stored for reference)
        y               % y coordinate matrix (stored for reference)
    end

    methods
        function obj = HermiteBeam(x, y, params)
            % Constructor
            % x, y:   Cartesian coordinate matrices
            % params: HermiteParameters object

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
            % Reconstructs HermiteParameters at the requested z to correctly
            % evaluate waist w(z) and Gouy phase psi(z) before computing the field.
            w0       = obj.Parameters.InitialWaist;
            lambda   = obj.Lambda;
            n        = obj.Parameters.n;
            m        = obj.Parameters.m;
            params_z = HermiteParameters(z, w0, lambda, n, m);
            field    = obj.computeField(X, Y, params_z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            w0     = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            params = GaussianParameters(z, w0, lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'hermite_3_2'.
            name = sprintf('hermite_%d_%d', obj.Parameters.n, obj.Parameters.m);
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, x, y, params)
            % Compute HG field from Cartesian grids and parameters.
            r  = sqrt(x.^2 + y.^2);
            GB = GaussianBeam(r, params);

            w  = params.Waist;
            n  = params.n;
            m  = params.m;

            % Hermite polynomial scaling: sqrt(2)*coord/w
            Hn = PolynomialUtils.hermitePoly(n, sqrt(2) * x ./ w);
            Hm = PolynomialUtils.hermitePoly(m, sqrt(2) * y ./ w);

            % PhiPhase = (n+m)*psi(z) — modal Gouy phase shift
            field = Hn .* Hm .* exp(1i * params.PhiPhase) .* GB.OpticalField;
        end
    end
end
