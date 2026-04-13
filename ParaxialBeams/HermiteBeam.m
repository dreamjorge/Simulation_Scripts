classdef HermiteBeam < ParaxialBeam
    % HermiteBeam - Scalar optical field for a Hermite-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = HermiteBeam(w0, lambda, n, m)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % e.g. 'hermite_3_2'
    %
    % Mathematical definition (standard Hermite-Gaussian, HG_{nm}):
    %
    %   u_{nm}(x,y,z) = H_n(sqrt(2)*x/w) * H_m(sqrt(2)*y/w)
    %                   * u_0(r,z) * exp(i*(n+m)*psi(z))
    %
    % where:
    %   H_n, H_m  - Hermite polynomials of order n, m
    %   w = w(z)  - beam waist at z
    %   u_0(r,z)  - fundamental Gaussian carrier field
    %   psi(z)    - Gouy phase = arctan(z / z_R)
    %
    % The coordinate system is Cartesian (x, y). The class stores only the
    % physical parameters w0, lambda, n, m — no grid or stored field.
    %
    % Reference: Saleh & Teich, "Fundamentals of Photonics", Ch. 3.

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        n               % Hermite order in x
        m               % Hermite order in y
    end

    methods
        function obj = HermiteBeam(w0, lambda, n, m)
            % Constructor
            % w0:     initial beam waist at z = 0 (m)
            % lambda: wavelength (m)
            % n, m:   Hermite mode orders (default 0)

            if nargin > 0
                obj = obj@ParaxialBeam(lambda);
                obj.InitialWaist = w0;
                if nargin >= 4
                    obj.n = n;
                    obj.m = m;
                else
                    obj.n = 0;
                    obj.m = 0;
                end
            else
                obj = obj@ParaxialBeam();
                obj.n = 0;
                obj.m = 0;
            end
        end

        % -----------------------------------------------------------------
        % ParaxialBeam interface
        % -----------------------------------------------------------------

        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex optical field on Cartesian grid (X,Y) at depth z.
            field = obj.computeField(X, Y, z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            params = GaussianParameters(z, obj.InitialWaist, obj.Lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'hermite_3_2'.
            name = sprintf('hermite_%d_%d', obj.n, obj.m);
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, X, Y, z)
            % Compute HG_{nm} field at Cartesian grid (X,Y) and depth z.
            w0     = obj.InitialWaist;
            lambda = obj.Lambda;
            k      = obj.k;
            zr     = pi * w0^2 / lambda;

            w   = w0 * sqrt(1 + (z/zr)^2);
            Rc  = z  * (1 + (zr/z)^2);
            if z == 0, Rc = Inf; end
            psi = atan2(z, zr);

            % Gaussian carrier field u_0(r,z)
            r          = sqrt(X.^2 + Y.^2);
            amplitude  = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
            phase_z    = -1i * k * z;
            phase_curv = 1i * k * r.^2 ./ (2 * Rc);
            phase_curv(isinf(Rc)) = 0;
            phase_gouy = -1i * psi;
            carrier    = amplitude .* exp(phase_z + phase_curv + phase_gouy);

            % Hermite polynomial scaling: sqrt(2)*coord/w
            Hn = PolynomialUtils.hermitePoly(obj.n, sqrt(2) * X ./ w);
            Hm = PolynomialUtils.hermitePoly(obj.m, sqrt(2) * Y ./ w);

            % Modal Gouy phase shift: (n+m)*psi
            phi_mode = (obj.n + obj.m) * psi;

            field = Hn .* Hm .* exp(1i * phi_mode) .* carrier;
        end
    end
end
