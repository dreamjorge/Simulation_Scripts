classdef GaussianBeam < ParaxialBeam
    % GaussianBeam - Scalar optical field for a fundamental Gaussian beam
    % Compatible with GNU Octave and MATLAB
    %
    % Constructor (Phase 3 API):
    %   beam = GaussianBeam(w0, lambda)
    %
    % Usage:
    %   field = beam.opticalField(X, Y, z)    % complex field on Cartesian grid
    %   params = beam.getParameters(z)         % GaussianParameters at z
    %   name   = beam.beamName()               % 'gaussian'
    %
    % The beam stores only the physical parameters w0 and lambda. All grid
    % quantities (field, waist at z, radius, Gouy phase) are computed on
    % demand by opticalField(X, Y, z). This makes the object stateless with
    % respect to propagation distance.

    properties
        InitialWaist    % Beam waist at z = 0 (m)
        OpticalField    % Legacy snapshot field compatibility
    end

    methods
        function obj = GaussianBeam(arg1, arg2, varargin)
            % Constructor
            % Modern API:
            %   GaussianBeam(w0, lambda)
            %
            % Legacy-compatible APIs:
            %   GaussianBeam(R, gaussianParams)
            %   GaussianBeam(X, Y, gaussianParams)

            if nargin == 0
                obj = obj@ParaxialBeam();
                obj.OpticalField = [];
                return;
            end

            if nargin == 2 && isa(arg2, 'GaussianParameters')
                params = arg2;
                obj = obj@ParaxialBeam(params.Lambda);
                obj.InitialWaist = params.InitialWaist;
                obj.OpticalField = obj.opticalFieldLegacy(arg1, [], params.zCoordinate);
                return;
            end

            if nargin == 3 && isa(varargin{1}, 'GaussianParameters')
                params = varargin{1};
                obj = obj@ParaxialBeam(params.Lambda);
                obj.InitialWaist = params.InitialWaist;
                obj.OpticalField = obj.opticalFieldLegacy(arg1, arg2, params.zCoordinate);
                return;
            end

            if nargin >= 2
                w0 = arg1;
                lambda = arg2;
                obj = obj@ParaxialBeam(lambda);
                obj.InitialWaist = w0;
                obj.OpticalField = [];
            else
                obj = obj@ParaxialBeam();
                obj.OpticalField = [];
            end
        end

        % -----------------------------------------------------------------
        % ParaxialBeam interface
        % -----------------------------------------------------------------

        function field = opticalField(obj, X, Y, z)
            % opticalField - Complex optical field on Cartesian grid (X,Y) at depth z.
            %
            % Field formula (fundamental Gaussian TEM_{00}):
            %
            %   u(r,z) = (w0/w) * exp(-r^2/w^2)
            %            * exp(-i*k*z) * exp(i*k*r^2/(2R)) * exp(-i*psi)
            %
            % where:
            %   r   = sqrt(x^2 + y^2)    radial distance (m)
            %   w   = w(z)                beam waist at z
            %   R   = R(z)                radius of curvature (Inf at z=0)
            %   psi = psi(z) = arctan(z/zR)   Gouy phase

            w0     = obj.InitialWaist;
            lambda = obj.Lambda;
            k      = obj.k;
            zr     = pi * w0^2 / lambda;

            w   = w0 * sqrt(1 + (z ./ zr).^2);
            Rc  = z  .* (1 + (zr ./ z).^2);
            if isscalar(z) && z == 0
                Rc = Inf;
            else
                Rc(z == 0) = Inf;
            end
            psi = atan2(z, zr);

            r          = sqrt(X.^2 + Y.^2);
            amplitude  = (w0 ./ w) .* exp(-r.^2 ./ w.^2);
            phase_z    = -1i * k .* z;
            phase_curv = 1i * k .* r.^2 ./ (2 .* Rc);
            phase_curv(isinf(Rc)) = 0;
            phase_gouy = -1i * psi;

            field = amplitude .* exp(phase_z + phase_curv + phase_gouy);
        end

        function params = getParameters(obj, z)
            % getParameters - GaussianParameters evaluated at axial position z.
            params = GaussianParameters(z, obj.InitialWaist, obj.Lambda);
        end

        function name = beamName(obj)
            name = 'gaussian';
        end
    end

    methods (Access = private)
        function field = opticalFieldLegacy(obj, coordA, coordB, z)
            if isempty(coordB)
                X = coordA;
                Y = zeros(size(coordA));
            else
                X = coordA;
                Y = coordB;
            end
            field = obj.opticalField(X, Y, z);
        end
    end
end
