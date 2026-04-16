classdef LaguerreBeam < ParaxialBeam
    % LaguerreBeam - Scalar optical field for a Laguerre-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = LaguerreBeam(w0, lambda, l, p)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % e.g. 'laguerre_2_1'
    %
    % Mathematical definition (standard Laguerre-Gaussian, LG_{lp}):
    %
    %   u_{lp}(r,theta,z) = (sqrt(2)*r/w)^|l| * L_p^|l|(2*r^2/w^2)
    %                       * u_0(r,z) * exp(i*l*theta) * exp(i*(|l|+2p)*psi(z))
    %
    % where:
    %   l         - topological charge (azimuthal index, integer, may be negative)
    %   p         - radial order (non-negative integer)
    %   L_p^|l|   - associated Laguerre polynomial of degree p and order |l|
    %   w = w(z)  - beam waist at z
    %   u_0(r,z)  - fundamental Gaussian carrier field
    %   psi(z)    - Gouy phase = arctan(z / z_R)
    %
    % The class accepts Cartesian (X, Y) in opticalField and converts to polar
    % internally. The class stores only w0, lambda, l, p — no grid or field.
    %
    % Reference: Allen et al., PRA 45, 8185 (1992).

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        l               % Topological charge (azimuthal index)
        p               % Radial order
        OpticalField    % Legacy snapshot field compatibility
    end

    methods
        function obj = LaguerreBeam(arg1, arg2, varargin)
            % Constructor
            % Modern API:
            %   LaguerreBeam(w0, lambda, l, p)
            %
            % Legacy-compatible API:
            %   LaguerreBeam(R, Theta, laguerreParams)

            % Call superclass constructor first (MATLAB requirement)
            obj = obj@ParaxialBeam();

            % Determine parameters from input using static helper
            [w0, lambda, l, p, legacyCoords, legacyZ] = ...
                LaguerreBeam.parseArgs(arg1, arg2, varargin{:});

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
            % Converts Cartesian to polar internally, then evaluates the LG field.
            [TH, R] = cart2pol(X, Y);
            field   = obj.computeField(R, TH, z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            params = GaussianParameters(z, obj.InitialWaist, obj.Lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'laguerre_2_1'.
            name = sprintf('laguerre_%d_%d', obj.l, obj.p);
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

            if nargin == 3 && isa(varargin{1}, 'LaguerreParameters')
                % Legacy: LaguerreBeam(R, Theta, laguerreParams)
                params = varargin{1};
                lambda = params.Lambda;
                w0 = params.InitialWaist;
                l = params.l;
                p = params.p;
                legacyCoords{1} = arg1;
                legacyCoords{2} = arg2;
                legacyZ = params.zCoordinate;

            elseif nargin >= 2
                % Modern: LaguerreBeam(w0, lambda, l, p)
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
            % Compute LG_{lp} field from polar grids (r, theta) and depth z.
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

            % Radial amplitude: (sqrt(2)*r/w)^|l|
            amp_lg = (sqrt(2) * r ./ w).^abs(l);

            % Associated Laguerre polynomial: L_p^|l|(2*r^2/w^2)
            Lpl = PolynomialUtils.associatedLaguerre(p, l, 2 * r.^2 ./ w.^2);

            % Modal Gouy phase shift: (|l|+2p)*psi
            phi_mode = (abs(l) + 2*p) * psi;

            field = amp_lg .* Lpl .* exp(1i * l * theta) ...
                    .* exp(-1i * phi_mode) .* carrier;
        end
    end
end
