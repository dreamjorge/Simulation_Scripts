classdef HankelLaguerre < ParaxialBeam
    % HankelLaguerre - Hankel-type Laguerre-Gaussian beam field
    % Compatible with GNU Octave and MATLAB
    %
    % Coordinate system: internally POLAR (r, theta).
    % The unified API opticalField(X, Y, z) accepts Cartesian and converts
    % internally via cart2pol.
    %
    % -------------------------------------------------------------------------
    % Theoretical Context (Contexto Academico)
    % -------------------------------------------------------------------------
    % Hankel beams are specialized solutions to the paraxial wave equation
    % that represent "conical" waves. Unlike standard LG modes, these beams
    % possess a net radial energy flow, making them ideal for modelling
    % diverging or converging wavefronts in complex media.
    %
    % Mathematical Definition:
    % They are constructed through the superposition of a paraxial mode
    % and its quadrature (Hilbert transform) companion:
    %
    %   H_{lp}^{(1,2)}(r, phi, z) = LG_{lp}(r, phi, z) +/- i * XLG_{lp}(r, phi, z)
    %
    % where XLG_{lp} is the Hilbert-transformed counterpart defined as:
    %   - Azimuthal phase: exp(-i*p*theta) instead of exp(i*l*theta)
    %   - Radial/polynomial structure: same as LG_{lp}
    %
    %   XLG_{lp}(r, phi, z) = amp(r,z) * L_p^|l|(x) * exp(-i*p*theta) *
    %                          exp(i*(|l|+2p)*psi(z)) * u_0(r,z)
    %
    % Physical Significance:
    %   H^{(1)}: Outward propagating wave (away from optical axis) — net radial
    %            energy flux pointing outward, analogous to Hankel function H^(1).
    %   H^{(2)}: Inward  propagating wave (towards optical axis) — net radial
    %            energy flux pointing inward, analogous to Hankel function H^(2).
    %
    % The construction is analogous to combining real and imaginary parts of
    % the scalar Hankel functions in cylindrical wave theory:
    %   H^{(1)}(kr) = J(kr) + i*Y(kr)
    %
    % Reference: Kotlyar et al., "Hankel-Bessel laser beams", J. Opt. Soc. Am. A.

    properties
        Parameters           % LaguerreParameters object (at construction z)
        HankelType           % 1 = H^(1) outward, 2 = H^(2) inward
        OpticalFieldLaguerre % Complex field array (kept for backward compatibility)
    end

    methods
        function obj = HankelLaguerre(r, theta, params, hankelType)
            % Constructor
            % r, theta:   polar coordinate matrices
            % params:     LaguerreParameters object
            % hankelType: 1 (H^(1), default) or 2 (H^(2))

            if nargin > 0
                if nargin < 4
                    hankelType = 1;
                end
                obj = obj@ParaxialBeam(params.Lambda);
                obj.Parameters           = params;
                obj.HankelType           = hankelType;
                obj.OpticalFieldLaguerre = hankelField(r, theta, params, hankelType);
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
            % Converts Cartesian to polar, reconstructs LaguerreParameters at z,
            % then computes the Hankel superposition.
            [TH, R]  = cart2pol(X, Y);
            w0       = obj.Parameters.InitialWaist;
            lambda   = obj.Lambda;
            l        = obj.Parameters.l;
            p        = obj.Parameters.p;
            params_z = LaguerreParameters(z, w0, lambda, l, p);
            field    = hankelField(R, TH, params_z, obj.HankelType);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            w0     = obj.Parameters.InitialWaist;
            lambda = obj.Lambda;
            params = GaussianParameters(z, w0, lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'hankel1_laguerre_2_1'.
            name = sprintf('hankel%d_laguerre_%d_%d', ...
                obj.HankelType, obj.Parameters.l, obj.Parameters.p);
        end
    end
end

% -------------------------------------------------------------------------
% Module-private helper (not a method — keeps the class definition clean)
% -------------------------------------------------------------------------
function field = hankelField(r, theta, params, hankelType)
    % hankelField - Compute H^(1) or H^(2) Hankel-Laguerre field.
    %
    % Both components (LG and XLG) share the same radial amplitude and
    % Laguerre polynomial; they differ only in the azimuthal phase:
    %   LG : exp( i*l*theta)
    %   XLG: exp(-i*p*theta)   <- Hilbert-transform companion

    l   = params.l;
    p   = params.p;
    w   = params.Waist;

    % Standard LG amplitude term: (sqrt(2)*r/w)^|l|
    amp  = (sqrt(2) * r ./ w).^abs(l);

    % Associated Laguerre polynomial: L_p^|l|(2*r^2/w^2)
    xArg = 2 * r.^2 ./ w.^2;
    Lpl  = PolynomialUtils.associatedLaguerre(p, l, xArg);

    % Gaussian carrier field (common to both components)
    GB     = GaussianBeam(r, params);
    GField = GB.OpticalField;

    % Standard LG field
    LB_field  = amp .* Lpl .* exp( 1i * l * theta) .* exp(1i * params.PhiPhase) .* GField;

    % XLG (Hilbert-transformed) field — same amplitude, phase uses -p instead of l
    XLG_field = amp .* Lpl .* exp(-1i * p * theta) .* exp(1i * params.PhiPhase) .* GField;

    % Hankel superposition
    if hankelType == 1
        field = LB_field + 1i * XLG_field;   % H^(1): outward
    else
        field = LB_field - 1i * XLG_field;   % H^(2): inward
    end
end
