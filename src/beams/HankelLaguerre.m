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
        OpticalFieldLaguerre % Legacy compatibility field container
    end

    methods
        function obj = HankelLaguerre(w0, lambda, l, p, hankelType, varargin)
            % Constructor
            % w0:         initial beam waist at z = 0 (m)
            % lambda:     wavelength (m)
            % l:          topological charge (default 0)
            % p:          radial order (default 0)
            % hankelType: 1 = H^(1) [default], 2 = H^(2)
            %
            % Legacy constructor (also supported):
            %   HankelLaguerre(r, theta, laguerreParameters, hankelType)
            % Produces .OpticalFieldLaguerre for legacy scripts.
            %
            % Note: uses varargin to capture the raw 4th argument in legacy calls,
            % avoiding Octave parse errors when hankelType is not passed.

            % Call superclass constructor first (MATLAB requirement)
            obj = obj@ParaxialBeam();

            % Detect legacy API by checking the type of the 3rd argument.
            % Legacy: HankelLaguerre(r, theta, laguerreParams, hankelType)
            %   — arg3 is a LaguerreParameters object
            % Modern: HankelLaguerre(w0, lambda, l, p, hankelType)
            %   — arg3 is numeric (topological charge) or empty
            % Note: inputname(3) is empty when arg3 is a computed expression
            % (e.g. LaguerreParameters(...)), so we rely on isa() directly.
            if nargin >= 3 && isa(l, 'LaguerreParameters')
                isLegacy = true;
            else
                isLegacy = false;
            end

            if isLegacy
                % Legacy API: HankelLaguerre(r, theta, laguerreParams, hankelType)
                % The 4th legacy arg maps to 'p' in this signature (not varargin).
                if nargin >= 4
                    raw_hankelType_arg = p;
                else
                    raw_hankelType_arg = 1;
                end
                % parseArgs expects hankelType defined; set from extracted value.
                hankelType = raw_hankelType_arg;
            else
                % Modern API.
                raw_hankelType_arg = [];
                if nargin < 5, hankelType = 1; end
                if nargin < 4, p = 0; end
                if nargin < 3, l = 0; end
            end

            % Determine parameters from input using static helper
            [w0_out, lambda_out, l_out, p_out, hankelType_out, legacyCoords, legacyZ] = ...
                HankelLaguerre.parseArgs(w0, lambda, l, p, hankelType, raw_hankelType_arg);

            % Initialize parent class state
            if ~isempty(lambda_out)
                obj.Lambda = lambda_out;
                obj.k = 2 * pi / lambda_out;
            end

            % Initialize subclass state
            obj.InitialWaist = w0_out;
            obj.l = l_out;
            obj.p = p_out;
            obj.HankelType = hankelType_out;

            if ~isempty(legacyCoords{1})
                obj.OpticalFieldLaguerre = hankelField(legacyCoords{1}, legacyCoords{2}, ...
                    w0_out, lambda_out, l_out, p_out, legacyZ, hankelType_out);
            else
                obj.OpticalFieldLaguerre = [];
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

    methods (Static)
        function [w0, lambda, l, p, hankelType, legacyCoords, legacyZ] = parseArgs(w0, lambda, l, p, hankelType, raw_hankelType_arg)
            % Static helper to parse constructor arguments
            % Detects legacy API: when l is a LaguerreParameters object
            if nargin < 6, raw_hankelType_arg = []; end
            legacyCoords = {[], []};
            legacyZ = 0;

            if nargin < 1
                l = 0; p = 0; hankelType = 1;
                return;
            end

            % Legacy path: HankelLaguerre(r, theta, laguerreParams, hankelType)
            if nargin >= 3 && isa(l, 'LaguerreParameters')
                laguerreParams = l;
                legacyCoords{1} = w0;      % r coordinate
                legacyCoords{2} = lambda;  % theta coordinate
                legacyZ = laguerreParams.zCoordinate;
                w0 = laguerreParams.InitialWaist;
                lambda = laguerreParams.Lambda;
                l = laguerreParams.l;
                p = laguerreParams.p;
                if ~isempty(raw_hankelType_arg)
                    hankelType = raw_hankelType_arg;
                else
                    hankelType = 1;
                end
            else
                % Modern path: HankelLaguerre(w0, lambda, l, p, hankelType)
                if nargin < 3, l = 0; end
                if nargin < 4, p = 0; end
                if nargin < 5, hankelType = 1; end
            end
        end

        function Rays = getPropagateCylindricalRays(Rays, TotalRays, r, th, difr, LParametersZi, LParametersZ, HankelType)
            % getPropagateCylindricalRays - Legacy-compatible ray propagation API.
            tempdr = num2cell(difr);
            [dr, dtheta, dz] = deal(tempdr{:});

            for point_index = 1:TotalRays
                ri = Rays.rCoordinate(point_index) + (1 ./ Rays.zrSlope(point_index)) * dz;
                thi = Rays.thetaCoordinate(point_index) + (1 ./ Rays.zthSlope(point_index)) * dz;
                zi = Rays.zCoordinate(point_index) + dz;

                rayHankelType = HankelType;
                if (ri < 0) && (HankelType == 2)
                    rayHankelType = 1;
                    ri = abs(ri);
                end

                Rays.rCoordinate(point_index) = ri;
                Rays.thetaCoordinate(point_index) = thi;
                Rays.zCoordinate(point_index) = zi;
                Rays.hankelType(point_index) = rayHankelType;

                [xi, yi] = pol2cart(thi, ri);
                if (HankelType == 2) && (rayHankelType == 1)
                    xi = -xi;
                    yi = -yi;
                end
                Rays.xCoordinate(point_index) = xi;
                Rays.yCoordinate(point_index) = yi;

                HLr = HankelLaguerre(r, thi, LParametersZi, rayHankelType);
                HLth = HankelLaguerre(ri, th, LParametersZi, rayHankelType);
                HLz = HankelLaguerre(ri, thi, LParametersZ, rayHankelType);

                fr = unwrap(angle(HLr.OpticalFieldLaguerre));
                fth = unwrap(angle(HLth.OpticalFieldLaguerre));
                fz = unwrap(angle(HLz.OpticalFieldLaguerre));

                [zrSlope, zthSlope, rthSlope] = HankelLaguerre.gradientCylindrical(fr, fth, fz, LParametersZi.k, dr, dtheta, dz, ri, thi, zi);
                Rays.zrSlope(point_index) = zrSlope;
                Rays.zthSlope(point_index) = zthSlope;
                Rays.rthSlope(point_index) = rthSlope;
            end
        end
    end

    methods (Static, Access = private)
        function [zrSlope, zthSlope, rthSlope] = gradientCylindrical(fr, ftheta, fz, k, dr, dtheta, dz, r, theta, z)
            gr = gradient(fr) ./ dr;
            gtheta = (1 ./ r) .* (gradient(ftheta) ./ dtheta);
            gz = gradient(fz) ./ dz + k;

            n = size(gr, 2);
            idxZ = HankelLaguerre.clampIndex(floor(z / dz + 1), numel(gz));
            idxR = HankelLaguerre.clampIndex(n / 2 + 1 + floor(r / dr), numel(gr));
            idxTheta = HankelLaguerre.clampIndex(n / 2 + 1 + floor(theta / dtheta), numel(gtheta));

            zrSlope = gz(idxZ) / gr(idxR);
            zthSlope = gz(idxZ) / gtheta(idxTheta);
            rthSlope = gr(idxR) / gtheta(idxTheta);
        end

        function idx = clampIndex(raw, maxSize)
            idx = floor(raw);
            if idx < 1
                idx = 1;
            elseif idx > maxSize
                idx = maxSize;
            end
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
    LB_field  = amp .* Lpl .* exp( 1i * l * theta) .* exp(-1i * phi_mode) .* carrier;

    % XLG (Hilbert-transformed) field — same amplitude, phase uses -p instead of l
    XLG_field = amp .* Lpl .* exp(-1i * p * theta) .* exp(-1i * phi_mode) .* carrier;

    % Hankel superposition
    if hankelType == 1
        field = LB_field + 1i * XLG_field;   % H^(1): outward
    else
        field = LB_field - 1i * XLG_field;   % H^(2): inward
    end
end
