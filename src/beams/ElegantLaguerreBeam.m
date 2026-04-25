classdef ElegantLaguerreBeam < ParaxialBeam
    % ElegantLaguerreBeam - Elegant Laguerre-Gaussian beam implementation
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = ElegantLaguerreBeam(w0, lambda, l, p)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % e.g. 'elegant_laguerre_2_1'
    %
    % Mathematical differences from standard LaguerreBeam:
    %
    %   Standard LG: amplitude = (sqrt(2)*r/w)^|l|,  poly_arg = 2*r^2/w^2
    %   Elegant LG:  amplitude = (sqrt(alpha)*r)^|l|, poly_arg = alpha*r^2
    %
    % where alpha(z) = i*k / (2*q(z)), q(z) = z + i*z_R (complex beam parameter).
    %
    % Full field definition (elegant Laguerre-Gaussian, ELG_{lp}):
    %
    %   u_{lp}(r,theta,z) = (sqrt(alpha)*r)^|l| * L_p^|l|(alpha*r^2)
    %                       * u_0(r,z) * exp(i*l*theta) * exp(i*(|l|+2p)*psi(z))
    %
    % Physical consequence: using the complex beam parameter alpha instead of
    % the real waist w changes both the radial amplitude profile and the
    % polynomial argument, producing the "elegant" amplitude-phase coupling.
    %
    % The class accepts Cartesian (X, Y) in opticalField and converts to polar
    % internally. Stores only w0, lambda, l, p — no grid or field.
    %
    % Reference: Siegman, "Lasers", University Science Books (1986), Ch. 17;
    %            Siegman, J. Opt. Soc. Am. A 13, 952 (1996).

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        l               % Topological charge (azimuthal index)
        p               % Radial order
        OpticalField    % Legacy snapshot field compatibility
    end

    methods
        function obj = ElegantLaguerreBeam(arg1, arg2, varargin)
            % Constructor
            % Modern API:
            %   ElegantLaguerreBeam(w0, lambda, l, p)
            %
            % Legacy-compatible API:
            %   ElegantLaguerreBeam(R, Theta, laguerreParams)

            % Call superclass constructor first (MATLAB requirement)
            obj = obj@ParaxialBeam();

            % Emit deprecation warning (Strangler Fig migration)
            warning('BeamFactory:deprecated', ...
                'src/beams/ElegantLaguerreBeam is deprecated. Use BeamFactory.create(''elegant_laguerre'', ...) or +paraxial/+beams/ElegantLaguerreBeam directly.');

            % Determine parameters from input using static helper
            [w0, lambda, l, p, legacyCoords, legacyZ] = ...
                ElegantLaguerreBeam.parseArgs(arg1, arg2, varargin{:});

            % Initialize parent class state
            if ~isempty(lambda)
                obj.Lambda = lambda;
                obj.k = 2 * pi / lambda;
            end

            % Initialize subclass state
            obj.InitialWaist = w0;
            obj.l = l;
            obj.p = p;

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
            %
            % Converts Cartesian to polar internally, then evaluates the ELG field.
            [TH, R] = cart2pol(X, Y);
            field   = obj.computeField(R, TH, z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            params = GaussianParameters(z, obj.InitialWaist, obj.Lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'elegant_laguerre_2_1'.
            name = sprintf('elegant_laguerre_%d_%d', obj.l, obj.p);
        end
    end

    methods (Static)
        function [w0, lambda, l, p, legacyCoords, legacyZ] = parseArgs(arg1, arg2, varargin)
            % Static helper to parse constructor arguments
            w0 = [];
            lambda = [];
            l = 0;
            p = 0;
            legacyCoords = {[], []};
            legacyZ = 0;

            if nargin < 2
                return;
            end

            if nargin == 3 && (isa(varargin{1}, 'ElegantLaguerreParameters') || isa(varargin{1}, 'LaguerreParameters'))
                % Legacy: ElegantLaguerreBeam(R, Theta, laguerreParams)
                params = varargin{1};
                lambda = params.Lambda;
                w0 = params.InitialWaist;
                l = params.l;
                p = params.p;
                legacyCoords{1} = arg1;
                legacyCoords{2} = arg2;
                legacyZ = params.zCoordinate;

            elseif nargin >= 2
                % Modern: ElegantLaguerreBeam(w0, lambda, l, p)
                w0 = arg1;
                lambda = arg2;
                if numel(varargin) >= 2
                    l = varargin{1};
                    p = varargin{2};
                end
            end
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, r, theta, z)
            % Compute ELG_{lp} field from polar grids (r, theta) and depth z.
            w0     = obj.InitialWaist;
            lambda = obj.Lambda;
            k      = obj.k;
            l      = obj.l;
            p      = obj.p;
            zr     = pi * w0^2 / lambda;

            w   = w0 * sqrt(1 + (z/zr)^2);
            Rc  = z  * (1 + (zr/z)^2);
            if z == 0, Rc = Inf; end
            psi = atan2(z, zr);

            % Gaussian carrier field u_0(r,z)
            amplitude  = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
            phase_z    = -1i * k * z;
            phase_curv = 1i * k * r.^2 ./ (2 * Rc);
            phase_curv(isinf(Rc)) = 0;
            phase_gouy = -1i * psi;
            carrier    = amplitude .* exp(phase_z + phase_curv + phase_gouy);

            % Complex beam parameter: alpha(z) = i*k / (2*q(z))
            q         = z + 1i * zr;
            alpha_val = 1i * k / (2 * q);

            % Polynomial argument: alpha*r^2  (complex, unlike 2*r^2/w^2 in standard LG)
            Lpl = PolynomialUtils.associatedLaguerre(p, l, alpha_val .* r.^2);

            % Radial amplitude: (sqrt(alpha)*r)^|l|  -- note: alpha is complex
            amp_el = (sqrt(alpha_val) .* r).^abs(l);

            % Modal Gouy phase shift: (|l|+2p)*psi
            phi_mode = (abs(l) + 2*p) * psi;

            field = amp_el .* Lpl .* exp(1i * l * theta) ...
                    .* exp(-1i * phi_mode) .* carrier;
        end
    end
end
