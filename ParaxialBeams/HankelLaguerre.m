classdef HankelLaguerre < ParaxialBeam
    % HankelLaguerre - Hankel-type Laguerre-Gaussian beam field
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = HankelLaguerre(w0, lambda, l, p, hankelType)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % e.g. 'hankel1_laguerre_2_1'
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
    % where XLG_{lp} is the Hilbert-transformed counterpart:
    %   - Azimuthal phase: exp(-i*p*theta) instead of exp(i*l*theta)
    %   - Radial/polynomial structure: same as LG_{lp}
    %
    % Physical Significance:
    %   H^{(1)}: Outward propagating wave (away from axis) — analogous to H^(1)(kr)
    %   H^{(2)}: Inward  propagating wave (towards axis)  — analogous to H^(2)(kr)
    %
    % Reference: Kotlyar et al., "Hankel-Bessel laser beams", J. Opt. Soc. Am. A.

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        l               % Topological charge (azimuthal index)
        p               % Radial order
        HankelType      % 1 = H^(1) outward, 2 = H^(2) inward
    end

    methods
        function obj = HankelLaguerre(w0, lambda, l, p, hankelType)
            % Constructor
            % w0:         initial beam waist at z = 0 (m)
            % lambda:     wavelength (m)
            % l:          topological charge (default 0)
            % p:          radial order (default 0)
            % hankelType: 1 = H^(1) [default], 2 = H^(2)

            if nargin > 0
                obj = obj@ParaxialBeam(lambda);
                obj.InitialWaist = w0;
                if nargin >= 4
                    obj.l = l;
                    obj.p = p;
                else
                    obj.l = 0;
                    obj.p = 0;
                end
                if nargin >= 5
                    obj.HankelType = hankelType;
                else
                    obj.HankelType = 1;
                end
            else
                obj = obj@ParaxialBeam();
                obj.l          = 0;
                obj.p          = 0;
                obj.HankelType = 1;
            end
        end

        % -----------------------------------------------------------------
        % ParaxialBeam interface
        % -----------------------------------------------------------------

        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex optical field on Cartesian grid (X,Y) at depth z.
            %
            % Converts Cartesian to polar, then computes the Hankel superposition.
            [TH, R] = cart2pol(X, Y);
            field   = hankelField(R, TH, obj.InitialWaist, obj.Lambda, obj.l, obj.p, z, obj.HankelType);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            params = GaussianParameters(z, obj.InitialWaist, obj.Lambda);
        end

        function name = beamName(obj)
            % beamName - Returns identifier string, e.g. 'hankel1_laguerre_2_1'.
            name = sprintf('hankel%d_laguerre_%d_%d', obj.HankelType, obj.l, obj.p);
        end
    end
end

% -------------------------------------------------------------------------
% Module-private helper (not a method — keeps the class definition clean)
% -------------------------------------------------------------------------
function field = hankelField(r, theta, w0, lambda, l, p, z, hankelType)
    % hankelField - Compute H^(1) or H^(2) Hankel-Laguerre field at depth z.
    %
    % Both components (LG and XLG) share the same radial amplitude and
    % Laguerre polynomial; they differ only in the azimuthal phase:
    %   LG : exp( i*l*theta)
    %   XLG: exp(-i*p*theta)   <- Hilbert-transform companion

    k  = 2*pi/lambda;
    zr = pi * w0^2 / lambda;

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

    % Standard LG amplitude: (sqrt(2)*r/w)^|l|
    amp  = (sqrt(2) * r ./ w).^abs(l);

    % Associated Laguerre polynomial: L_p^|l|(2*r^2/w^2)
    Lpl  = PolynomialUtils.associatedLaguerre(p, l, 2 * r.^2 ./ w.^2);

    % Modal Gouy phase shift: (|l|+2p)*psi
    phi_mode = (abs(l) + 2*p) * psi;

    % Standard LG field
    LB_field  = amp .* Lpl .* exp( 1i * l * theta) .* exp(1i * phi_mode) .* carrier;

    % XLG (Hilbert-transformed) field — same amplitude, phase uses -p instead of l
    XLG_field = amp .* Lpl .* exp(-1i * p * theta) .* exp(1i * phi_mode) .* carrier;

    % Hankel superposition
    if hankelType == 1
        field = LB_field + 1i * XLG_field;   % H^(1): outward
    else
        field = LB_field - 1i * XLG_field;   % H^(2): inward
    end
end
