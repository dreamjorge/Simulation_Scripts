classdef ElegantLaguerreBeam < ParaxialBeam
    % ElegantLaguerreBeam - Elegant Laguerre-Gaussian beam implementation
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: internally POLAR (r, theta) - same as LaguerreBeam.
    % The unified API opticalField(X, Y, z) accepts Cartesian and converts
    % internally via cart2pol.
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
    %   u_{lp}(r,theta,z) = (sqrt(alpha)*r)^|l| * L_p^|l|(alpha*r^2) *
    %                       u_0(r,z) * exp(i*l*theta) * exp(i*(|l|+2p)*psi(z))
    %
    % Physical consequence: using the complex beam parameter alpha instead of
    % the real waist w changes both the radial amplitude profile and the
    % polynomial argument, producing the "elegant" amplitude-phase coupling.
    %
    % Reference: Siegman, "Lasers", University Science Books (1986), Ch. 17;
    %            Siegman, J. Opt. Soc. Am. A 13, 952 (1996).

    properties
        Parameters      % ElegantLaguerreParameters object
        OpticalField    % 2D array representing the field (at constructor z)
        r               % radial coordinate matrix (stored for reference)
        theta           % angular coordinate matrix (stored for reference)
    end

    methods
        function obj = ElegantLaguerreBeam(r, theta, params)
            % Constructor
            % r, theta: polar coordinate matrices
            % params:   ElegantLaguerreParameters object

            if nargin > 0
                obj = obj@ParaxialBeam(params.Lambda);
                obj.Parameters   = params;
                obj.r            = r;
                obj.theta        = theta;
                obj.OpticalField = obj.computeField(r, theta, params);
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
            % Converts Cartesian to polar, reconstructs ElegantLaguerreParameters
            % at z (recomputes alpha(z)), then evaluates the ELG field.
            [TH, R]  = cart2pol(X, Y);
            w0       = obj.Parameters.InitialWaist;
            lambda   = obj.Lambda;
            l        = obj.Parameters.l;
            p        = obj.Parameters.p;
            params_z = ElegantLaguerreParameters(z, w0, lambda, l, p);
            field    = obj.computeField(R, TH, params_z);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            w0     = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            params = GaussianParameters(z, w0, lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'elegant_laguerre_2_1'.
            name = sprintf('elegant_laguerre_%d_%d', obj.Parameters.l, obj.Parameters.p);
        end
    end

    % -----------------------------------------------------------------
    % Private helpers
    % -----------------------------------------------------------------
    methods (Access = private)
        function field = computeField(obj, r, theta, params)
            % Compute ELG field from polar grids and parameters.
            alpha_val = params.alpha;

            % Polynomial argument: alpha*r^2 (complex, unlike 2*r^2/w^2 in standard LG)
            xArg = alpha_val .* r.^2;
            Lpl  = PolynomialUtils.associatedLaguerre(params.p, params.l, xArg);

            GB = GaussianBeam(r, params);

            % Radial amplitude: (sqrt(alpha)*r)^|l|  — note: alpha is complex
            amplitudeTerm = (sqrt(alpha_val) .* r).^abs(params.l);

            % PhiPhase = (|l|+2p)*psi(z) — same Gouy shift as standard LG
            field = amplitudeTerm .* Lpl .* exp(1i * params.l * theta) ...
                    .* exp(1i * params.PhiPhase) .* GB.OpticalField;
        end
    end
end
