classdef LaguerreBeam < ParaxialBeam
    % LaguerreBeam - Scalar optical field for a Laguerre-Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: internally POLAR (r, theta).
    % The unified API opticalField(X, Y, z) accepts Cartesian coordinates
    % and converts to polar internally via cart2pol.
    %
    % Mathematical definition (standard Laguerre-Gaussian, LG_{lp}):
    %
    %   u_{lp}(r,theta,z) = (sqrt(2)*r/w)^|l| * L_p^|l|(2*r^2/w^2) *
    %                       u_0(r,z) * exp(i*l*theta) * exp(i*(|l|+2p)*psi(z))
    %
    % where:
    %   l         - topological charge (azimuthal index, integer, can be negative)
    %   p         - radial order (non-negative integer)
    %   L_p^|l|   - associated Laguerre polynomial of degree p and order |l|
    %   w = w(z)  - beam waist at z (from LaguerreParameters)
    %   u_0(r,z)  - fundamental Gaussian carrier field
    %   psi(z)    - Gouy phase = arctan(z / z_R)
    %
    % The azimuthal phase exp(i*l*theta) produces an optical vortex of charge l.
    % The extra Gouy phase (|l|+2p)*psi is in LaguerreParameters.PhiPhase.
    %
    % Reference: Allen et al., PRA 45, 8185 (1992).

    properties
        Parameters      % LaguerreParameters object
        OpticalField    % 2D array representing the field (at constructor z)
        r               % radial coordinate matrix (stored for reference)
        theta           % angular coordinate matrix (stored for reference)
    end

    methods
        function obj = LaguerreBeam(r, theta, params)
            % Constructor
            % r, theta: polar coordinate matrices
            % params:   LaguerreParameters object

            obj = obj@ParaxialBeam(params.Lambda);

            if nargin > 0
                obj.Parameters   = params;
                obj.r            = r;
                obj.theta        = theta;
                obj.OpticalField = obj.computeField(r, theta, params);
            end
        end

        % -----------------------------------------------------------------
        % ParaxialBeam interface
        % -----------------------------------------------------------------

        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex optical field on Cartesian grid (X,Y) at depth z.
            %
            % Converts Cartesian to polar, reconstructs LaguerreParameters at z,
            % then evaluates the LG field.
            [TH, R]  = cart2pol(X, Y);
            w0       = obj.Parameters.InitialWaist;
            lambda   = obj.Lambda;
            l        = obj.Parameters.l;
            p        = obj.Parameters.p;
            params_z = LaguerreParameters(z, w0, lambda, l, p);
            field    = obj.computeField(R, TH, params_z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            w0     = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            params = GaussianParameters(z, w0, lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'laguerre_2_1'.
            name = sprintf('laguerre_%d_%d', obj.Parameters.l, obj.Parameters.p);
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, r, theta, params)
            % Compute LG field from polar grids and parameters.
            GB = GaussianBeam(r, params);

            w   = params.Waist;
            l   = params.l;
            p   = params.p;

            % Radial amplitude: (sqrt(2)*r/w)^|l|
            amplitudeTerm = (sqrt(2) * r ./ w).^abs(l);

            % Associated Laguerre polynomial: L_p^|l|(2*r^2/w^2)
            xArg = 2 * r.^2 ./ w.^2;
            Lpl  = PolynomialUtils.associatedLaguerre(p, l, xArg);

            % PhiPhase = (|l|+2p)*psi(z) — modal Gouy phase shift
            field = amplitudeTerm .* Lpl .* exp(1i * l * theta) ...
                    .* exp(1i * params.PhiPhase) .* GB.OpticalField;
        end
    end
end
