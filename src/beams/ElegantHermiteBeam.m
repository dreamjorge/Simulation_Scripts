classdef ElegantHermiteBeam < ParaxialBeam
    % ElegantHermiteBeam - Elegant Hermite-Gaussian beam implementation
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = ElegantHermiteBeam(w0, lambda, n, m)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % e.g. 'elegant_hermite_2_0'
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
    %   u_{nm}(x,y,z) = H_n(sqrt(alpha)*x) * H_m(sqrt(alpha)*y)
    %                   * u_0(r,z) * exp(i*(n+m)*psi(z))
    %
    % Physical consequence: because alpha is complex, the Hermite polynomials
    % evaluated at complex arguments produce amplitude AND phase modulation
    % simultaneously — the hallmark of the "elegant" variant.
    %
    % Reference: Siegman, "Lasers", University Science Books (1986), Ch. 17;
    %            Siegman, J. Opt. Soc. Am. A 13, 952 (1996).

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        n               % Horizontal mode index
        m               % Vertical mode index
        OpticalField    % Legacy snapshot field compatibility
    end

    methods
        function obj = ElegantHermiteBeam(arg1, arg2, varargin)
            % Constructor
            % Modern API:
            %   ElegantHermiteBeam(w0, lambda, n, m)
            %
            % Legacy-compatible API:
            %   ElegantHermiteBeam(X, Y, hermiteParams)

            % Call superclass constructor first (MATLAB requirement)
            obj = obj@ParaxialBeam();

            % Emit deprecation warning (Strangler Fig migration)
            warning('BeamFactory:deprecated', ...
                'src/beams/ElegantHermiteBeam is deprecated. Use BeamFactory.create(''elegant_hermite'', ...) or +paraxial/+beams/ElegantHermiteBeam directly.');

            % Determine parameters from input using static helper
            [w0, lambda, n, m, legacyCoords, legacyZ] = ...
                ElegantHermiteBeam.parseArgs(arg1, arg2, varargin{:});

            % Initialize parent class state
            if ~isempty(lambda)
                obj.Lambda = lambda;
                obj.k = 2 * pi / lambda;
            end

            % Initialize subclass state
            obj.InitialWaist = w0;
            obj.n = n;
            obj.m = m;

            if ~isempty(legacyCoords{1})
                obj.OpticalField = obj.computeField(legacyCoords{1}, legacyCoords{2}, legacyZ);
            else
                obj.OpticalField = [];
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
            % beamName - Returns identifier string, e.g. 'elegant_hermite_2_0'.
            name = sprintf('elegant_hermite_%d_%d', obj.n, obj.m);
        end
    end

    methods (Static)
        function [w0, lambda, n, m, legacyCoords, legacyZ] = parseArgs(arg1, arg2, varargin)
            % Static helper to parse constructor arguments
            w0 = [];
            lambda = [];
            n = 0;
            m = 0;
            legacyCoords = {[], []};
            legacyZ = 0;

            if nargin < 2
                return;
            end

            if nargin == 3 && (isa(varargin{1}, 'ElegantHermiteParameters') || isa(varargin{1}, 'HermiteParameters'))
                % Legacy: ElegantHermiteBeam(X, Y, hermiteParams)
                params = varargin{1};
                lambda = params.Lambda;
                w0 = params.InitialWaist;
                n = params.n;
                m = params.m;
                legacyCoords{1} = arg1;
                legacyCoords{2} = arg2;
                legacyZ = params.zCoordinate;

            elseif nargin >= 2
                % Modern: ElegantHermiteBeam(w0, lambda, n, m)
                w0 = arg1;
                lambda = arg2;
                if numel(varargin) >= 2
                    n = varargin{1};
                    m = varargin{2};
                end
            end
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, X, Y, z)
            % Compute EHG_{nm} field at Cartesian grid (X,Y) and depth z.
            %
            % Key: sqrt(alpha) is complex, so H_n(sqrt(alpha)*x) is evaluated
            % with a complex argument, encoding both amplitude and phase.
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

            % Complex beam parameter: alpha(z) = i*k / (2*q(z))
            q          = z + 1i * zr;
            alpha      = 1i * k / (2 * q);
            sqrt_alpha = sqrt(alpha);

            % Hermite polynomials with complex argument
            Hn = PolynomialUtils.hermitePoly(obj.n, sqrt_alpha .* X);
            Hm = PolynomialUtils.hermitePoly(obj.m, sqrt_alpha .* Y);

            % Modal Gouy phase shift: (n+m)*psi
            phi_mode = (obj.n + obj.m) * psi;

            field = Hn .* Hm .* exp(-1i * phi_mode) .* carrier;
        end
    end
end
